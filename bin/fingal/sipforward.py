from esys.escript import *
from esys.escript.pdetools import ArithmeticTuple, PCG
from .tools import setupERTPDE

class SIPSolver(object):
    """
    a solver complex electrical conductivity problems using Schur complement.
    """
    def __init__(self, dom, surface_mask = None,
                 rtol=1e-8, atol=0, pdetol=1e-10, iter_max=100, verbose=False):
        self.dom = dom
        self.rtol = rtol
        self.atol = atol
        self.iter_max = iter_max
        self.verbose = verbose
        self.surface_mask = surface_mask
        self.alpha = None
        # forward PDE:
        self.pde_fw = setupERTPDE(dom, tolerance=pdetol, isComplex=False)
        # preconditioner PDE
        self.pde_prec = setupERTPDE(dom, tolerance=pdetol, isComplex=False)

    def setSource(self, source_location= [0., 0., 0.], u_src=None, tag=None, sigma_src=1):
        """
        this sets the location of the source and source potential, that is the
        electric potential due to a point source at source_location, injection current I=1
        and real valued (!) conductivity sigma_src. if u_src is not given the respective PDE is solved.
        Note: this function also sets the surface factor alpha which needs the source location
        """
        self.source_location = source_location
        self.sigma_src = sigma_src
        # surface factor
        n=self.dom.getNormal()
        r = n.getX() - self.source_location
        self.alpha = inner(n, r) / length(r) ** 2
        if self.surface_mask:
            self.alpha *= (1 - self.surface_mask)
        self.pde_fw.setValue(A=self.sigma_src * kronecker(3), d=self.sigma_src * self.alpha)
        S = Scalar(0., DiracDeltaFunctions(self.dom))
        S.setTaggedValue(tag, 1)
        self.pde_fw.setValue(y_dirac=S, X = Data(), y=Data())
        self.u_src = self.pde_fw.getSolution()
        self.pde_fw.setValue(A=Data(), d=Data(), y_dirac=Data())
        if self.verbose:
            print("Source location  : ", self.source_location)
            print("       potential : ", self.u_src)
    def setPrimaryPotential(self, sigma_p, current = 1 ):
        """
        sets the primary potential using the source potential

        u_p = u_src/sigma_p *sigma_src * current

        """
        self.u_p = self.u_src/sigma_p * self.sigma_src * current
        self.sigma_p = sigma_p
        self.sigma_bc_p = sigma_p
        self.current = current
        if self.verbose:
            print("primary conductivity  : ", self.sigma_p)
            print("        potential real : ", self.u_p.real())
            print("                  imag : ", self.u_p.imag())

    def setConductivity(self, sigma, sigma_bc = None):
        """
        set the conductivity sigma (interior) and sigma_bc (boundary).
        """
        sigma_re = sigma.real()
        assert inf(sigma_re) > 0, "real part of sigma must be positive."
        self.sigma_re = sigma_re
        self.sigma_im = sigma.imag()

        if sigma_bc is None:
            sigma_bc = interpolate(sigma, self.alpha.getFunctionSpace())
        else:
            sigma_bc = interpolate(sigma_bc, self.alpha.getFunctionSpace())
        sigma_bc_re = sigma_bc.real()
        assert inf(sigma_bc_re) > 0, "real part of sigma must be positive on the boundary."
        self.sigma_bc_re = sigma_bc_re
        self.sigma_bc_im = sigma_bc.imag()

        self.sigma_prec = self.sigma_re * ( 1 + (self.sigma_im / self.sigma_re) ** 2 )
        self.sigma_bc_prec = self.sigma_bc_re * ( 1 + (self.sigma_bc_im / self.sigma_bc_re) ** 2 )

        self.pde_fw.setValue(A=self.sigma_re * kronecker(3), d=sigma_bc_re * self.alpha)
        self.pde_prec.setValue(A=self.sigma_prec * kronecker(3), d=self.sigma_bc_prec * self.alpha)
        if self.verbose:
            print("conductivity real : ", self.sigma_re)
            print("             imag : ", self.sigma_im)
            print("conductivity BC real : ", self.sigma_bc_re)
            print("                imag : ", self.sigma_bc_im)
            print("preconditioner    : ", self.sigma_prec)
            print("               BC : ", self.sigma_bc_prec)
        return self

    def _evalA(self,p_re):
        p_bc_re = interpolate(p_re, self.alpha.getFunctionSpace())
        self.pde_fw.setValue(X = - self.sigma_im * grad(p_re),
                             y =  - self.alpha * self.sigma_bc_im * p_bc_re,
                             y_dirac = Data() )
        p_im = self.pde_fw.getSolution()
        p_bc_im = interpolate(p_im, self.alpha.getFunctionSpace())
        A = self.sigma_re * grad(p_re) - self.sigma_im * grad(p_im)
        a = self.alpha * (self.sigma_bc_re * p_bc_re - self.sigma_bc_im * p_bc_im)
        return ArithmeticTuple(A, a)

    def _evalM(self, R):
        self.pde_prec.setValue(X = R[0], y=R[1],  y_dirac = Data() )
        p_re = self.pde_prec.getSolution()
        return p_re

    def _bilinearform(self,  p_re,  R):
        return integrate(inner(R[0], grad(p_re))) + integrate(inner(R[1], p_re))

    def getU_im(self, u_re):

        self.pde_fw.setValue(X=(self.sigma_p * grad(self.u_p)).imag()  - self.sigma_im * grad(u_re),
                             y=self.alpha * ( (self.sigma_bc_p * self.u_p).imag() - self.sigma_bc_im * u_re),
                             y_dirac = Data())
        u_im = self.pde_fw.getSolution()
        return u_im

    def solve(self, u0_re=None):
        if u0_re is None:
            u0_re = self.u_p.real()
        u0_im=self.getU_im(u0_re)
        #u_re = u_inital.real()
        #u_im = u_inital.real()
        R0 = (self.sigma_p * grad(self.u_p)).real() - self.sigma_re * grad(u0_re) + self.sigma_im * grad(u0_im)
        r0 = self.alpha * ( (self.sigma_bc_p * self.u_p).real() - self.sigma_bc_re * u0_re + self.sigma_bc_im * u0_im)
        self.R_0 = ArithmeticTuple(R0, r0)

        u_re, R_k, history = PCG(self.R_0, self._evalA, u0_re, self._evalM, self._bilinearform, atol=self.atol,
                rtol=self.rtol, iter_max=self.iter_max, initial_guess=True, return_history = True, verbose=self.verbose)
        self.R_k = R_k
        self.history = history
        u_im = self.getU_im(u_re)
        return u_re + 1j * u_im