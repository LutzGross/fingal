from esys.escript import *
from fingal import setupERTPDE
import numpy as np
from numpy.linalg import norm


class MineGeometry(object):
    """
    class to old some infos
    """
    Stations = {}
    ElectrodeSpacing = None
    NumElectrodes = None
    LineOffset = None
    ExtractionWidth = None
    LineHeight = None
    TagFaces = None
    RemainderLength = None
    RoadHeight = None
    LineNorth = []
    LineSouth = []

def getGeometryFromGeoFile(geofile):
    out=MineGeometry()
    for line in open(geofile, 'r'):
        if line.strip().startswith("NumElectrodes"):
            out.NumElectrodes = int(line.split('=')[1][:-2])
        elif line.strip().startswith("LineOffset"):
            out.LineOffset = int(line.split('=')[1][:-2])
        elif line.strip().startswith("ElectrodeSpacing"):
            out.ElectrodeSpacing = float(line.split('=')[1][:-2])
        elif line.strip().startswith("ExtractionWidth"):
            out.ExtractionWidth = float(line.split('=')[1][:-2])
        elif line.strip().startswith("LineHeight"):
            out.LineHeight = float(line.split('=')[1][:-2])
        elif line.strip().startswith("RemainderLength"):
            out.RemainderLength = float(line.split('=')[1][:-2])
        elif line.strip().startswith("RoadHeight"):
            out.RoadHeight = float(line.split('=')[1][:-2])
        elif line.strip().startswith("Physical Surface"):
            out.TagFaces = line.split('"')[1]
    print(f"Dimension of electrode array read from file {geofile}.")
    print(f"NumElectrodes= {out.NumElectrodes}.")
    print(f"SpacingElectrodes = {out.ElectrodeSpacing}.")
    print(f"LineOffset = {out.LineOffset}.")
    print(f"LineHeight = {out.LineHeight}.")
    print(f"ExtractionWidth = {out.ExtractionWidth}.")
    print(f"RemainderLength = {out.RemainderLength}.")
    print(f"RoadHeight = {out.RoadHeight}.")
    print(f"Tag for faces = '{out.TagFaces}'.")

    # south line:
    for k in range(out.NumElectrodes):
        stid = 100 + (k+1)
        out.LineSouth.append(stid)
        out.Stations[stid]=np.array([out.LineOffset + k * out.ElectrodeSpacing, -out.ExtractionWidth / 2, out.LineHeight])
    # north line:
    for k in range(out.NumElectrodes):
        stid = 200 + (k+1)
        out.LineNorth.append(stid)
        out.Stations[stid]=np.array([out.LineOffset + k * out.ElectrodeSpacing, out.ExtractionWidth / 2, out.LineHeight])
    print(f"{len(out.Stations)} electodes found.")
    return out

def makePrimaryPotentials(domain, minegeo, sigma_ref, survey):
    injections=[]
    for A, B in survey.keys():
        if not A in injections:
            injections.append(A)
        if not B in injections:
            injections.append(B)
    mask_faces = Scalar(0, FunctionOnBoundary(domain))
    mask_faces.setTaggedValue(minegeo.TagFaces, 1)
    n = domain.getNormal() * mask_faces
    x = n.getX()
    primary_potentials = {}
    pde = setupERTPDE(domain, tolerance=1e-10)
    pde.setValue(A=sigma_ref * kronecker(3), y_dirac=Data(), X=Data(), Y=Data(), y=Data())
    for A in injections:
        s = Scalar(0., DiracDeltaFunctions(domain))
        s.setTaggedValue(f"s{A}", 1.)
        pde.setValue(y_dirac=s)
        r = x - minegeo.Stations[A]
        pde.setValue(d=sigma_ref * inner(r, n) / length(r) ** 2)  # doi:10.1190/1.1440975
        primary_potentials[A] = pde.getSolution()
        assert Lsup(primary_potentials[A]) > 0, "Zero potential for injection %s" % A
        print(f"primary potential for injection at {A}: {str(primary_potentials[A])}.")
    print(f"{len(primary_potentials)} primary potentials calculated.")
    return primary_potentials

def makeSecondaryPotentials(domain, minegeo, sigma, sigma_ref, primary_potentials):
    """
    It is assumed that sigma_ref=sigma on faces!!!!
    """

    mask_faces = Scalar(0, FunctionOnBoundary(domain))
    mask_faces.setTaggedValue(minegeo.TagFaces, 1)
    n = domain.getNormal() * mask_faces
    x = n.getX()
    pde = setupERTPDE(domain, tolerance=1e-10)
    secondary_potentials = {}
    pde.setValue(A=sigma * kronecker(3), y_dirac=Data(), X=Data(), Y=Data(), y=Data())
    for A in primary_potentials:
        r = x - minegeo.Stations[A]
        pde.setValue(d=sigma_ref * inner(r, n) / length(r) ** 2)
        pde.setValue(X=(sigma_ref-sigma)*grad(primary_potentials[A]))
        secondary_potentials[A] = pde.getSolution()
        print(f"secondary potential for injection at {A}: {str(secondary_potentials[A])}.")
    print(f"{len(secondary_potentials)} secondary potentials calculated.")
    return secondary_potentials


def makeMeasurements(line, data, minegeo, injections=[], dir=0):
    M = []
    X = []
    for k in range(len(line)-1):
        if not ( line[k] in injections or line[k+1] in injections ):
            m=data[k]-data[k+1]
            x=(minegeo.Stations[line[k]][dir]+minegeo.Stations[line[k+1]][dir])/2
            X.append(x)
            M.append(m)
    return X, M

def makeResistivity1(domain, width_fracture_zone, rho_ref, rho_frac, minegeo):
    global rho
    # X=ReducedFunction(domain).getX()
    X = Function(domain).getX()
    # fractures are west of excavated zone
    kk = 1
    m1 = whereNonPositive(X[0] - minegeo.RemainderLength)
    m2 = clip((X[0] - (minegeo.RemainderLength - width_fracture_zone - minegeo.FringFractureZone)) / minegeo.FringFractureZone,
              minval=0, maxval=1)
    d1 = clip((X[1] - (-minegeo.ExtractionWidth / 2 + minegeo.ResetFactureZoneSouth)) / minegeo.FringFractureZone, minval=0, maxval=1)
    d2 = clip(((minegeo.ExtractionWidth / 2 - minegeo.ResetFactureZoneNorth) - X[1]) / minegeo.FringFractureZone, minval=0, maxval=1)
    e1 = clip((X[2] - (minegeo.RoadHeight - minegeo.FringFractureZone)) / minegeo.FringFractureZone, minval=0, maxval=1)
    e2 = clip(((minegeo.RoadHeight + minegeo.ThicknessFactureZone + minegeo.FringFractureZone) - X[2]) / minegeo.FringFractureZone, minval=0,
              maxval=1)
    f = m1 * m2 * d1 * d2 * e1 * e2
    rho = rho_ref * (1 - f) + rho_frac * f
    return rho

def makeApparentResitivity(line, data, minegeo, injections=(), I =1, dir=0):
    A, B = injections
    XA=minegeo.Stations[A]
    XB=minegeo.Stations[B]
    RHO = []
    X = []
    for k in range(len(line)-1):
        if not ( line[k] in injections or line[k+1] in injections ):
            XM = minegeo.Stations[line[k]]
            XN = minegeo.Stations[line[k+1]]
            g = 1/(4*np.pi) * (1/norm(XM-XA)-1/norm(XM-XB)-1/norm(XN-XA)+1/norm(XN-XB))
            rho=(data[k]-data[k+1])/(I*g)
            x=(minegeo.Stations[line[k]][dir]+minegeo.Stations[line[k+1]][dir])/2
            X.append(x)
            RHO.append(rho)

    return X, RHO