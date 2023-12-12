#!/usr/bin/python3
from esys.escript import *
import importlib, sys, os
sys.path.append(os.getcwd())
import numpy as np
from .datahandling import SurveyData
from fingal import readElectrodeLocations, readSurveyData, setupERTPDE, getSourcePotentials, makeMaskForOuterSurface

from esys.escript.pdetools import Locator, MaskFromTag



class IPSynthetic(object):
    """
    This creates synthetic data files for ERT and IP for know conductivity and normalized chargeability
    distributions sigma_0 and M_n.  
    """
    def __init__(self, domain, schedule, sigma_src=1,  mask_faces=None, stationsFMT="e%s",
                 createSecondaryData=True, createFieldData=False, printInfo = True):
        """
        :param domain: physical domain
        :type domain: `AbstractDomain`
        :param schedule: schedule of the survey 'ABMN' (or AMN, etc)
        :type schedule: `SurveyData`
        :param mask_faces: mask of the faces with electric radiation condition. If None, the front, back, left, right and
                            bottom faces are used.
        :type mask_faces: `Data` of `FunctionOnBoundary` or None
        :param createSecondaryData: runs the synthetic survey for secondary voltages needed for
                                    chargeability data. set to false if these data are not created to save compute time
        :param createFieldData: runs the synthetic survey for fulL waver field data set to false if these data
                                are not created to save compute time
        :param sigma_src: conductivity value to create source potential to suppress singularities at sources.
        :param stationsFMT: format string to convert stations keys to mesh tags.
        :type schedule: `SurveyData`
        """
        assert isinstance(schedule, SurveyData)
        self.domain = domain
        self.printinfo = printInfo
        self.createFieldData = createFieldData
        self.createSecondaryData =createSecondaryData
        self.schedule = schedule
        self.sigma_src=sigma_src
        self.stationsFMT=stationsFMT
        self.mask_faces=mask_faces

        # 
        station_locations = [ schedule.getStationLocationByKey(S) for S in schedule.getStationNumeration() ]

        self.nodelocators = Locator(Solution(domain), station_locations)
        self.elementlocators = Locator(ReducedFunction(domain), station_locations)
        self.stationlocators = Locator(Function(domain), station_locations)

        self.source_potential = getSourcePotentials(domain, sigma_src, self.schedule, mask_outer_faces=mask_faces,
                                                    stationsFMT=self.stationsFMT)

        self.source_field = {}
        self.source_potential_at_station = {}
        self.source_field_at_station = {}
        self.sigma_src = {}

        for iA in self.source_potential:
            self.source_potential_at_station[iA] = self.nodelocators(self.source_potential[iA])
            self.sigma_src[iA] = sigma_src
            self.source_field[iA] = -grad(self.source_potential[iA], ReducedFunction(self.domain))
            self.source_field_at_station[iA] = self.elementlocators(self.source_field[iA])
            if self.printinfo:
                print("\tSource: %s : sigma = %s " % (iA, str(self.sigma_src[iA])))


        if self.printinfo:
            print("%s electrode locations found" % (len(self.sigma_src)))
            print(str(len(station_locations)) + " station locators calculated.")
            print(str(len(self.source_field)) + " source fields calculated.")

    def setProperties(self, sigma_0=1., sigma_0_faces=None, M_n=None, M_n_faces=None):
        """
        sets the DC conductivity sigma_0 and normalized chargeability.
        if createSecondaryData M_n and M_n_faces must be given.
        """

        if self.printinfo:
            print(".... DC potential : (increment from source)")
            print("sigma_0 = ", str(sigma_0) )

        # -div(sigma_src/alpha_A grad(U_A*alpha_A))=S_A
        # U_S = injection potential
        # secondary potential -div(sigma_0 grad(W_A) = -div( (sigma_src/alpha_A - sigma_0) grad(U_A*alpha_A))
        #  DC potenential V_0 = W_A + U_A*alpha_A)
        pde = setupERTPDE(self.domain)

        self.alpha={}
        self.potential_0 = {}
        self.potential_0_at_stations = {}
        if self.createFieldData:
            self.field_0_at_stations = {}
            self.field_0 = {}

        sigma_0_at_stations=self.stationlocators(sigma_0)
        n = self.domain.getNormal() * self.mask_faces
        x_bc = FunctionOnBoundary(self.domain).getX()

        pde.setValue(A=sigma_0 * kronecker(3), y_dirac=Data())
        for iA in self.source_potential:
            alpha_A = self.sigma_src[iA]/sigma_0_at_stations[iA]
            self.source_potential[iA]*= alpha_A
            self.sigma_src[iA]*= 1./alpha_A

            xA = self.schedule.getStationLocationByNumber(iA)
            r = x_bc - xA
            fA=inner(r, n) / length(r) ** 2
            pde.setValue(d=sigma_0_faces * fA, y= ( self.sigma_src[iA] - sigma_0_faces ) * fA * self.source_potential[iA])
            pde.setValue(X=(self.sigma_src[iA] - sigma_0) * grad(self.source_potential[iA]))
            self.potential_0[iA] = pde.getSolution()
            self.potential_0_at_stations[iA] = self.nodelocators(self.potential_0[iA])
            if self.createFieldData:
                self.field_0[iA] = -grad(self.potential_0[A], ReducedFunction(self.domain))
                self.field_0_at_stations[iA] = self.elementlocators(self.field_0[iA])
            if self.printinfo:
                print("\t%s : dV0 = %s " % (iA, str(self.potential_0[iA])))

        if self.printinfo:
            print(str(len(self.potential_0)) + " DC potentials calculated.")
            if self.createFieldData:
                print(str(len(self.field_0)) + " DC fields calculated.")

        if self.createSecondaryData:
            # secondary potential -div(sigma_oo grad(V_2) = -div(-Mn grad(V_0))
            if M_n is None or M_n_faces is None:
                raise ValueError("Secondary potential needed but no normalized chargeability M_n or M_n_faces give")
            sigma_oo = M_n + sigma_0
            sigma_oo_faces = M_n_faces + sigma_0_faces
            if self.printinfo:
                print(".... secondary (IP) potentials:")
                print("sigma_oo = ", str(sigma_oo))
                print("sigma_oo_faces = ", str(sigma_oo_faces))
                print("M_n = ", str(M_n))
                print("M_n_faces = ", str(M_n_faces))

            self.secondary_potential = {}
            self.secondary_potential_at_stations = {}
            if self.createFieldData:
                self.secondary_field = {}
                self.secondary_field_at_stations = {}
            pde.setValue(A=sigma_oo * kronecker(3), y_dirac=Data())
            for iA in self.source_potential:
                V_0=self.potential_0[iA] + self.source_potential[iA]
                xA = self.schedule.getStationLocationByNumber(iA)
                r = x_bc - xA
                fA = inner(r, n) / length(r) ** 2
                pde.setValue(d=sigma_oo_faces * fA, y=-M_n_faces * fA * V_0)
                pde.setValue(X=-M_n * grad(V_0))

                self.secondary_potential[iA] = pde.getSolution()
                self.secondary_potential_at_stations[iA] = self.nodelocators(self.secondary_potential[iA])
                if self.createFieldData:
                    self.secondary_field[iA] = -grad(self.secondary_potential[iA], ReducedFunction(self.domain))
                    self.secondary_field_at_stations[iA] = self.elementlocators(self.secondary_field[iA])
                if self.printinfo:
                    print("\t%s : V2 = %s " % (iA, str(self.secondary_potential[iA])))
            if self.printinfo:
                print(str(len(self.secondary_potential)) + " secondary potentials calculated.")
                if self.createFieldData:
                    print(str(len(self.secondary_field)) + " secondary fields calculated.")
    def write(self, filename, datacolumns = ['R'], addNoise = False,
                            rel_error=5, delimiter=",", usesStationCoordinates= False,
                            iFMT="%d", dFMT="%.5g", xFMT="%e"):
        """
        writes a data file with given datacolumns

        :param filename: name of data file to be created
        :param datacolumns: list of the data values to be created, see `SurveyData`
        :param noise: random noise in % to be added to the data
        :param rel_error: relative error in % used to write error columns when requested.
        """
        if rel_error > 0 :
            rel_error = rel_error / 100.
        else:
            rel_error = 0.05
        if self.printinfo:
            print(f"Data columns to be generated: {datacolumns}")
            if addNoise:
                print(f"relative error added to data = {rel_error}.")
            else:
                print(f"assumed relative error {rel_error} (not added).")

        dV_s = self.schedule.makeResistencePrediction(values=self.source_potential_at_station, valuesKeysAreStationKeys = False)
        dV_0 = self.schedule.makeResistencePrediction(values=self.potential_0_at_stations, valuesKeysAreStationKeys = False)
        dV_2 = self.schedule.makeResistencePrediction(values=self.secondary_potential_at_stations, valuesKeysAreStationKeys = False)
        if self.createFieldData:
            dE_s = self.schedule.makeResistencePrediction(values=self.source_field_at_station, valuesKeysAreStationKeys = False)
            dE_0 = self.schedule.makeResistencePrediction(values=self.field_0_at_stations, valuesKeysAreStationKeys = False)
            dE_2 = self.schedule.makeResistencePrediction(values=self.secondary_field_at_stations, valuesKeysAreStationKeys = False)

        #  now the data file can be created:
        if getMPIRankWorld() == 0:
            n = 0
            f = open(filename, 'w')

            FMTX = xFMT + delimiter + xFMT + delimiter + xFMT + delimiter
            FMTG = dFMT + delimiter + " "
            FMTI = iFMT + delimiter + " "
            for t in self.schedule.tokenIterator():

                out = ""
                if self.schedule.hasDipoleInjections():
                    A = t[0]
                    B = t[1]
                    if usesStationCoordinates:
                        out = FMTX % self.schedule.getStationLocationByKey(A) + FMTX % self.schedule.getStationLocationByKey(B)
                    else:
                        out = FMTI % A + FMTI % B
                    m = t[2:]
                    s = t[:2]
                else:
                    A = t[0]
                    s = t[:1]
                    m = t[1:]
                    if  usesStationCoordinates:
                        out = FMTX % self.schedule.getStationLocationByKey(A)
                    else:
                        out = FMTI % A
                if self.schedule.hasDipoleMeasurements():
                    M = m[0]
                    N = m[1]
                    if usesStationCoordinates:
                        out += FMTX % self.schedule.getStationLocationByKey(M) + FMTX % self.schedule.getStationLocationByKey(N)
                    else:
                        out += FMTI % M + FMTI % N
                else:
                    M = m[0]
                    if usesStationCoordinates:
                        out += FMTX % self.schedule.getStationLocationByKey(M)
                    else:
                        out += FMTI % M
                if addNoise:
                    pert = np.random.uniform(low = -rel_error, high = rel_error)
                else:
                    pert = 0
                V_0   = (dV_0[t] + dV_s[t]) * (1 + pert)  # ERT potential
                V_2   = dV_2[t] * (1 + pert)              # secondary potential
                V_oo  = V_0 + V_2
                ETA   = V_2 / V_0
                if self.createFieldData:
                    raise ValueError("CHECK E FIELD DATA")
                    E_0 = ( dE_0[t] + dE_s[t]) * (1 + pert)
                    E_2 = dE_2[t] * (1 + pert) # secondary potential
                    E_oo = E_0 + E_2
                    EI = sqrt(E_0[0] ** 2 + E_0[1] ** 2 + E_0[2] ** 2)
                    GAMMA = inner(E_2, E_0) / EI ** 2

                for o in datacolumns:
                    if o == 'R':  # resistance [V/A]
                        out += FMTG % (V_0)
                    elif o == 'E':  # electric field intensity per charge current [V/(Am)]
                        raise ValueError("check EI ")
                        out += FMTG % (EI)
                    elif o == 'ETA':  # chargeability  [1]
                        out += FMTG % (ETA )
                    elif o == 'GAMMA':  #
                        raise ValueError("check GAMMA ")
                        out += FMTG % (GAMMA )
                    elif o == 'ERR_R':  # resistance [V/A]
                        ERR_R = abs(V_0) * rel_error
                        out += FMTG % ERR_R
                    elif o == 'ERR_E':  # electric field intensity per charge current [V/(Am)]
                        raise ValueError("check ERR_E ")
                        ERR_E = abs(EI) * rel_error
                        out += FMTG % ERR_E
                    elif o == 'ERR_ETA':  # chargeability  [1]
                        ERR_ETA = 2 * abs(ETA) * rel_error
                        out += FMTG % ERR_ETA
                    elif o == 'ERR_GAMMA':  # modified chargeability (=ETA/(1-ETA)) [1]
                        raise ValueError("check gamma")
                        ERR_GAMMA = abs(GAMMA) * rel_error
                        out += FMTG % ERR_GAMMA
                    elif o == 'RELERR_R':  # resistance [1]
                        RELERR_R = rel_error
                        out += FMTG % RELERR_R
                    elif o == 'RELERR_E':  # electric field intensity per charge current [1]
                        raise ValueError("check RELERR_E ")
                        RELERR_E = rel_error
                        out += FMTG % RELERR_E
                    elif o == 'RELERR_ETA':  # chargeability  [1]
                        RELERR_ETA = 2 * rel_error
                        out += FMTG % RELERR_ETA
                    elif o == 'RELERR_GAMMA':  #
                        raise ValueError("check RELERR_GAMMA ")
                        RELERR_GAMMA = rel_error
                        out += FMTG % RELERR_GAMMA
                f.write(out[:-3] + "\n")
                n += 1
            f.close()
            if self.printinfo:
                print(n, " measurements were written to file " + filename)