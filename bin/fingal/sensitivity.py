from esys.escript import *
import importlib, sys, os
sys.path.append(os.getcwd())
import numpy as np
from .datahandling import SurveyData
from fingal import readElectrodeLocations, readSurveyData, setupERTPDE, getSourcePotentials, getAdditivePotentials

from esys.escript.pdetools import Locator, MaskFromTag



class ERTSensitivity(object):
    """
    This creates synthetic data files for ERT and IP for know conductivity and normalized chargeability
    distributions sigma_0 and M_n.  
    """
    def __init__(self, domain, schedule, sigma_src=1,  maskZeroPotential=None, stationsFMT="e%s",
                 createSecondaryData=True, createFieldData=False, printInfo = True):
        """
        :param domain: physical domain
        :type domain: `AbstractDomain`
        :param schedule: schedule of the survey 'ABMN' (or AMN, etc)
        :type schedule: `SurveyData`
        :param maskZeroPotential: mask where potentials are set to zero.
        :type maskZeroPotential: `Data` of `Solution` or None
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
        self.maskZeroPotential=maskZeroPotential

        # 
        station_locations = [ schedule.getStationLocationByKey(S) for S in schedule.getStationNumeration() ]

        self.nodelocators = Locator(Solution(domain), station_locations)
        self.elementlocators = Locator(ReducedFunction(domain), station_locations)
        self.stationlocators = Locator(Function(domain), station_locations)

        self.source_potential = getSourcePotentials(domain, sigma_src, self.schedule, maskZeroPotential=self.maskZeroPotential,
                                                    stationsFMT=self.stationsFMT)
        self.source_potential_at_station = { iA: self.nodelocators(self.source_potential[iA])  for iA in self.source_potential }

        self.potential_0_at_stations = None
        self.potential_0 = None

        if self.printinfo:
            print("source conductivity = %s"%str(self.sigma_src))
            print("%s electrode locations found" % (len(self.source_potential_at_station)))
            print(str(len(station_locations)) + " station locators calculated.")
            print(str(len(self.source_potential)) + " source potentials calculated.")


    def setProperties(self, sigma_0=1., sigma_0_at_stations=None):
        """
        sets the DC conductivity sigma_0
        """

        if self.printinfo:
            print(".... DC potential : (increment from source potential)")
            print("sigma_0 = ", str(sigma_0) )

        # -div(sigma_src/alpha_A grad(U_A*alpha_A))=S_A
        # U_S = injection potential
        # secondary potential -div(sigma_0 grad(V_0) = -div( (sigma_src/alpha_A - sigma_0) grad(U_A*alpha_A))
        #  DC potenential V_0 = W_A + U_A*alpha_A)
        pde = setupERTPDE(self.domain)
        pde.setValue(q=self.maskZeroPotential)
        if sigma_0_at_stations is None:
                sigma_0_at_stations = self.stationlocators(sigma_0)
        self.potential_0, src_potential_scale_0 = getAdditivePotentials(pde,
                                                  sigma = sigma_0,
                                                  schedule = self.schedule,
                                                  sigma_stations = sigma_0_at_stations,
                                                  source_potential = self.source_potential,
                                                  sigma_src = self.sigma_src)
        for iA in self.potential_0:
            self.potential_0[iA]+=self.source_potential[iA]*src_potential_scale_0[iA]

        self.potential_0_at_stations = { iA: self.nodelocators(self.potential_0[iA])  for iA in self.potential_0 }
        if self.printinfo:
            print(str(len(self.potential_0)) + " DC potentials calculated.")
            if self.createFieldData:
                print(str(len(self.field_0)) + " DC fields calculated.")
    def getSensitivity(self, *stations):
        if self.potential_0:
            potential_at_stations= self.potential_0_at_stations
            potential = self.potential_0
        else:
            potential_at_stations = self.source_potential_at_station
            potential = self.source_potential
        if self.schedule.hasDipoleInjections():
            if self.schedule.hasDipoleMeasurements():
                iA, iB, iM, iN = tuple( [ self.schedule.getStationNumber(ST) for ST in stations ] )
                F_ABMN = potential_at_stations[iA][iM]-potential_at_stations[iA][iN]-potential_at_stations[iB][iM]+potential_at_stations[iB][iN]
                F_ABMN_1 =  (self.source_potential_at_station[iA][iM]-self.source_potential_at_station[iA][iN]-self.source_potential_at_station[iB][iM]+self.source_potential_at_station[iB][iN]) * self.sigma_src
                UAB = potential[iA] - potential[iB]
                UMN = potential[iM] - potential[iN]
                sigma_a = F_ABMN_1/F_ABMN
                S_ABMN = sigma_a/F_ABMN * inner(grad(UAB, ReducedFunction(self.domain)), grad(UMN, ReducedFunction(self.domain)) )
                S_ABMN_max = Lsup(S_ABMN)
                print(stations, " -> sigma_a = ", sigma_a, "; S = ", S_ABMN_max )
            else:
                raise NotImplementedError()
        else:
            raise NotImplementedError()
        return S_ABMN

    def getTotalSensitivity(self):
        """
        get sensitivity across the survey.

        """
        out = None
        N=0
        for ST in self.schedule.tokenIterator():
            if N > 0:
                out+=abs(self.getSensitivity(*ST))
                N+=1
            else:
                out=abs(self.getSensitivity(*ST))
                N=1
        if N > 0:
            out*=1./N
        return out

