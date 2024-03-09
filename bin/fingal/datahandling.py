"""
fingal - data handling 

by l.gross@uq.edu.au, 2018
"""
import numpy as np
from esys.escript import getMPIRankWorld

def readElectrodeLocations(csvfile, delimiter=','):
    """
    
    reads the electrode locations from CSV  file in the form S, x, y, z
    and returns it as a dictionary { S : numpy.array((x,y,z))}
    """
    f=open(csvfile, 'r')
    locations={}
    line=f.readline().strip() 
    while line.strip():
        ll=line.strip().split(delimiter)
        S=int(ll[0])
        x=float(ll[1])
        y=float(ll[2])
        z=float(ll[3])
        locations[S]=np.array([x,y,z])
        line=f.readline()
    f.close()
    return locations    


def readSurveyData(csvfile, stations={}, usesStationCoordinates=False, columns=['R'], hasInjections=True, dipoleInjections=True, 
                           dipoleMeasurements=True, delimiter=',', commend='#', unDefined=-999999, printInfo=True):
    """
    creates a SurveyData object from a csvfile
    """
    nc=len(columns)
    if usesStationCoordinates:
        c=3
    else:
        c=1
    if hasInjections:
        if dipoleInjections and dipoleMeasurements:
            ns=4
        elif not dipoleInjections and dipoleMeasurements:
            ns=3
        elif dipoleInjections and not dipoleMeasurements:
            ns=3
        else:
            ns=2
    else:
        if dipoleMeasurements:
            ns=2
        else:
            ns=1
        
    f=open(csvfile, 'r')
    data=SurveyData(stations=stations, observations=columns, dipoleInjections=dipoleInjections, dipoleMeasurements=dipoleMeasurements, hasInjections=hasInjections, unDefined=unDefined)
    line=f.readline().strip()
    lc=0
    distmin=1e99
    distmax=0.
    while line:
        if not line.startswith(commend):
            ll=line.split(delimiter)
            if len(ll) < ns*c+nc:
                raise KeyError('In sufficient number of values in line %d. %d expected but %d found. line="%s"'%(lc, ns*c+nc, len(ll), line.strip()))
            if usesStationCoordinates:
                C=[ float(ll[i]) for i in range(ns*c) ]
                S=[]
                for i in range(ns): 
                    eid, edist = FindNearestElectrode(C[3*i], C[3*i+1], C[3*i+2], stations)
                    S.append(eid)
                    distmin=min(distmin, edist)
                    distmax=max(distmax, edist)
                #S=[ FindNearestElectrode(C[3*i], C[3*i+1], C[3*i+2], stations)[0] for i in range(ns) ]
            else:
                S=[ int(ll[i]) for i in range(ns) ]
            if ns==1:
                S=S[0]
            else:
                S=tuple(S)

            if nc >0:
                try:
                    F=[ float(ll[i+ns*c]) for i in range(nc) ]
                except ValueError:
                    raise ValueError(f"unable to read line {lc}: {line}")
                data.setDataRecord(S,tuple(F))
            else:
                data.setDataRecord(S,None)
            
        line=f.readline().strip()
        lc+=1
    f.close()
    if printInfo and getMPIRankWorld() == 0:
        print(f"Columns used: {columns}")
        print(f"{lc} data read from {csvfile}")
        if usesStationCoordinates:
            print("Maximum/Minimum distance of electrode locations to given stations = %e/%e."%(distmax, distmin))
    return data



class SurveyData(object):
    """
    this provides a simple interface to access GalvanicData type observations
    
        T->(columns[0], columns[1] ,...)
    
    where for injection electrodes A, B and observation electrodes M,N with the tokens 
            
            T=(A, B, M, N)
   
    or if dipoleInjections is not set 
            
            T=(A, B, M)

    or if dipoleMeasurements is not set 
            getSecondaryResistenceData
            T=(A, M, N)

    or if dipoleMeasurements and dipoleInjections are not set 
            
            T=(A, M)
    
    and (columns[0], columns[1], ....) refers to the data entries at each observation.
    
    
    The location of the injecting or observing electrode A is given by stationlocations[A]

    V_ABMN = V_A(X_M)-V_B(X_M)-V_A(X_N)-V_B(X_N)
    V_0 = voltage for sigma(omega=0)
    V_oo = voltage for sigma(omega=oo)
    V2 = V_oo-V_0 secondary voltage (or over-voltage)

    R = resistence V_0/I
    R2 = V2/I =  ETA * R = secondary resistance
    ETA = R2/R =V2/V_0 = intrinsic chargeabilityETA

    E = magnitude of electric field E0 of sigma(omega=0) divided by I
    E2 = (E_oo-E0).E0/|E0| = secondary electric field component (divided by I)
    """
    OBSTYPES= ['R',      'E',     'ETA',     'R2',     'E2',
               'ERR_R', 'ERR_E', 'ERR_ETA', 'ERR_R2', 'ERR_E2',
               'RELERR_R', 'RELERR_E', 'RELERR_ETA', 'RELERR_R2', 'RELERR_E2']
    def __init__(self, stations={}, observations=[], dipoleInjections=True, dipoleMeasurements=True,  hasInjections=True, default_rel_error=0.05, unDefined=-999999):
        """
        :stations: dictionary of station identifer to coordinates
        
        """
        for o in observations:
            assert self.checkObservationType(o), f'Unknown observation type {o}'
        self.default_rel_error=default_rel_error
        self.data={}
        self.stations=stations
        self.dipoleinjections=dipoleInjections and hasInjections
        self.dipolemeasurements=dipoleMeasurements
        self.hasinjections=hasInjections
        self.observations=observations
        self.lenObservations=len(observations)
        self.unDefined=unDefined
        if hasInjections:
            if dipoleInjections and dipoleMeasurements:
                ns=4
            elif not dipoleInjections and dipoleMeasurements:
                ns=3
            elif dipoleInjections and not dipoleMeasurements:
                ns=3
            else:
                ns=2
        else:
            if dipoleMeasurements:
                ns=2
            else:
                ns=1
        self.lenTokens=ns
        self.stations=stations
        
        # this makes sure that there is a fixed order being used:
        self.station_index={}
        self.station_key=[]
        for k in self.stations.keys():
            self.station_key.append(k)
            self.station_index[k]=self.station_key.index(k)
        self.injectionIter=None    # list of injection dipoles or stations
        self.injectionStations=None # list of injection stations
        self.observationElectrodes=None
        self._resistence_max=None
        self._secondary_resistence_max=None
    @classmethod 
    def checkObservationType(cls, obs):
        if isinstance(obs, list):
            return all([cls.checkObservationType(o) for o in obs])
        else:
            return obs in cls.OBSTYPES
    
    def isUndefined(self, v):
        if not v>self.unDefined+1:
            return True
        else:
            return False
        
    def getStationNumeration(self):
        """
        returns list of the station identifiers in a particular order
        """
        return self.station_key
    
    def getStationNumber(self, k):
        """
        returns the station number from station identifier key
        """
        return self.station_index[k]
    
    def getKeyOfStationNumber(self, ik):
        return self.station_key[ik]
    
    def hasDipoleInjections(self):
        return self.dipoleinjections
    
    def hasDipoleMeasurements(self):
        return self.dipolemeasurements
    def hasInjections(self):
        return self.hasinjections

    def getLengthTokens(self):
        return self.lenTokens
    def getLengthObservations(self):
        return self.lenObservations
    def getNumObservations(self):
        return len(self.data)
    def getNumStations(self):
        return len(self.station_key)
    def getStationLocation(self, k):
        return self.getStationLocationByKey(k)
    def getStationLocationByKey(self, k):
        return self.stations[k]
    def getStationLocationByNumber(self, ik):
        return self.getStationLocationByKey(self.getKeyOfStationNumber(ik))
    def setDataRecord(self, token, data=()):
        """
        sets observation by token.
        token needs to be a tuple of the form (A,B,M,N) (or simular) and data is a tuple giving the values of the observation
        """
        if self.getLengthTokens()>1:
            if not len(token) == self.getLengthTokens():
                raise ValueError('token length is incorrect.')
        if self.getLengthObservations() > 0:
            if not len(data) == self.getLengthObservations():
                raise ValueError('length of observations data is incorrect')
            self.data[token]=data
        else:
            self.data[token]=()
    def getDataRecord(self, token, datatype=None):
        """
        return the data for a given token (e.g. for (A,B,M,N) )
        """
        if datatype is None:
            return self.data[token]
        else:
            i=self.observations.index(datatype)
            return self.data[token][i]
    def hasDataType(self, datatype):
        return datatype in self.observations

    def getMaximumResistence(self):
        if self._resistence_max is None:
            self._resistence_max = max([ abs(self.getResistenceData(t)) for t in self.tokenIterator() ] )
        return self._resistence_max
    def getMaximumSecondaryResistence(self):
        if self._secondary_resistence_max is None:
            self._secondary_resistence_max = max([ abs(self.getSecondaryResistenceData(t)) for t in self.tokenIterator() ] )
        return self._secondary_resistence_max
    # ======================
    def getResistenceData(self, token):
        return self.getDataRecord(token, datatype='R')

    def getResistenceError(self, token):
        if self.hasDataType("ERR_R"):
            return self.getDataRecord(token, datatype='ERR_R')
        elif self.hasDataType("RELERR_R"):
            r=self.getDataRecord(token, datatype='R')
            e=self.getDataRecord(token, datatype='RELERR_R')
            if self.isUndefined(r) or self.isUndefined(e):
                return self.unDefined
            else:
                return abs(r*e)
        else:
            return self.default_rel_error * self.getMaximumResistence()
    def getResistenceRelError(self, token):
        if self.hasDataType("RELERR_R"):
            return self.getDataRecord(token, datatype='RELERR_R')
        elif self.hasDataType("ERR_R"):
            r=self.getDataRecord(token, datatype='R')
            e=self.getDataRecord(token, datatype='ERR_R')
            if self.isUndefined(r) or self.isUndefined(e):
                return self.unDefined
            elif abs(r)>0:
                return abs(e/r)
            else:    
                return self.unDefined
        else:
            return self.default_rel_error

    # ======================
    def getSecondaryResistenceData(self, token):
        if self.hasDataType("R2"):
            return self.getDataRecord(token, datatype='R2')
        else:
            eta=self.getChargeabilityData(token)
            r=self.getResistenceData(token)
            if self.isUndefined(r) or self.isUndefined(m):
                return self.unDefined
            else:
                return eta * r

    def getSecondaryResistenceError(self, token):
        if self.hasDataType("ERR_R2"):
            return self.getDataRecord(token, datatype='ERR_R2')
        elif self.hasDataType("ERR_ETA"):
            raise TypeError("REVISE")
            eta=self.getChargeabilityData(token)
            r=self.getResistenceData(token)
            eta_err=self.getChargeabilityError(token)
            r_err=self.getResistenceError(token)
            if self.isUndefined(r) or self.isUndefined(eta) or self.isUndefined(eta_err) or self.isUndefined(r_err) :
                return self.unDefined
            else:
                return r*eta_err+eta*r_err
        else:
            return self.default_rel_error * self.getMaximumSecondaryResistence()
    # =========================================
    def getFieldIntensityData(self, token):
        return self.getDataRecord(token, datatype='E')
    def getFieldIntensityError(self, token):
        if self.hasDataType("ERR_E"):
            return self.getDataRecord(token, datatype='ERR_E')
        elif self.hasDataType("RELERR_E"):
            r=self.getDataRecord(token, datatype='E')
            e=self.getDataRecord(token, datatype='RELERR_E')
            if self.isUndefined(r) or self.isUndefined(e):
                return self.unDefined
            else:
                return abs(e * r)
        else:
            return self.default_rel_error * self.getMaximumFieldIntensity()

    def getFieldIntensityRelError(self, token):
        if self.hasDataType("RELERR_E"):
            return self.getDataRecord(token, datatype='RELERR_E')
        elif self.hasDataType("ERR_E"):
            raise TypeError("REVISE")
            r=self.getDataRecord(token, datatype='E')
            e=self.getDataRecord(token, datatype='ERR_E')
            if self.isUndefined(r) or self.isUndefined(e):
                return self.unDefined
            elif abs(r)>0:
                return abs(e/r)
            else:    
                return self.unDefined
        else:
            return self.default_rel_error

    def tokenIterator(self):
        """
        returns a list of tokens (e.g. list of (A,B,M,N)s)
        """
        return self.data.keys()
    
    def injectionIterator(self):
        """
        returns an iterator (in form of list of id keys) of the injection dipols or stations:
        """
        if self.injectionIter is None:
            self.injectionIter=[]
            if self.hasInjections():
                if self.hasDipoleInjections() and self.hasDipoleMeasurements():
                    for A,B,M,N in self.tokenIterator():
                        if self.injectionIter.count((A,B)) ==0:
                            self.injectionIter.append((A,B))
                elif not self.hasDipoleInjections() and self.hasDipoleMeasurements():
                    for A,M,N in self.tokenIterator():
                        if self.injectionIter.count(A) ==0:
                            self.injectionIter.append(A)
                elif self.hasDipoleInjections() and not self.hasDipoleMeasurements():
                    for A,B, M in self.tokenIterator():
                        if self.injectionIter.count((A,B)) ==0:
                            self.injectionIter.append((A,B))
                else:
                    for A,M in self.tokenIterator():
                        if self.injectionIter.count(A) ==0:
                            self.injectionIter.append(A)
            else:
                self.injectionIter=[]
                
        return self.injectionIter
    def getInjectionStations(self):
        """
        returns a list of the stations used for injection
        """
        return self.getListOfInjectionStations()
    
    def getListOfInjectionStations(self):
        """
        returns a list of the stations keys used for injection
        """
        if self.injectionStations is None:
            out=[]
            if self.hasInjections():
                if self.hasDipoleInjections() and self.hasDipoleMeasurements():
                    for A,B,M,N in self.tokenIterator():
                        if out.count(B) ==0:
                            out.append(B)
                        if out.count(A) ==0:
                            out.append(A)

                elif not self.hasDipoleInjections() and self.hasDipoleMeasurements():
                    for A,M,N in self.tokenIterator():
                        if out.count(A) ==0:
                            out.append(A)
                elif self.hasDipoleInjections() and not self.hasDipoleMeasurements():
                    for A,B, M in self.tokenIterator():
                        if out.count(B) ==0:
                            out.append(B)
                        if out.count(A) ==0:
                            out.append(A)
                else:
                    for A,M in self.tokenIterator():
                        if out.count(A) ==0:
                            out.append(A)
            out.sort()
            self.injectionStations=out
        return self.injectionStations
    
    def getInjectionStationIndex(self, A):
        """
        returns the index of the injection station A in self.getListOfInjectionStations()
        """
        l=self.getListOfInjectionStations()
        if A in l:
            return l.index(A)
        else:
            raise ValueError("Unknown injection station %s"%A)

    def getListOfAllStations(self):
        """
        returns a list of the stations used for injection and observations (doubles are remove)
        """
        out=[]
        if self.hasInjections():
            if self.hasDipoleInjections() and self.hasDipoleMeasurements():
                for A,B,M,N in self.tokenIterator():
                    if out.count(B) ==0:
                        out.append(B)
                    if out.count(A) ==0:
                        out.append(A)
                    if out.count(M) ==0:
                        out.append(M)
                    if out.count(N) ==0:
                        out.append(N)

            elif not self.hasDipoleInjections() and self.hasDipoleMeasurements():
                for A,M,N in self.tokenIterator():
                    if out.count(A) ==0:
                        out.append(A)
                    if out.count(M) ==0:
                        out.append(M)
                    if out.count(N) ==0:
                        out.append(N)

            elif self.hasDipoleInjections() and not self.hasDipoleMeasurements():
                for A,B, M in self.tokenIterator():
                    if out.count(B) ==0:
                        out.append(B)
                    if out.count(A) ==0:
                        out.append(A)
                    if out.count(M) ==0:
                        out.append(M)
            else:
                for A,M in self.tokenIterator():
                    if out.count(A) ==0:
                        out.append(A)
                    if out.count(M) ==0:
                        out.append(M)

        return out          
    
    def getInjections(self, M=None, N=None):
        """
        returns a list of the injections that are measured by M,N 
        """
        out=[]
        if self.hasInjections():
            if self.hasDipoleInjections() and self.hasDipoleMeasurements():
                for A,B,M1,N1 in self.tokenIterator():
                    if M1 == M and N1 == N and out.count((A,B)) ==0:
                        out.append((A,B))

            elif not self.hasDipoleInjections() and self.hasDipoleMeasurements():
                for A,M1,N1 in self.tokenIterator():
                    if M1 == M and N1 == N and out.count(A) ==0:
                        out.append(A)
            elif self.hasDipoleInjections() and not self.hasDipoleMeasurements():
                for A,B, M1 in self.tokenIterator():
                    if M1 == M and out.count((A,B)) ==0:
                        out.append((A,B))
            else:
                for A,M1 in self.tokenIterator():
                    if M1 == M and out.count(A) ==0:
                        out.append(A)
        return out         

    def getObservationElectrodes(self):
        """
        returns a list of the observations electrodes:
        """
        if self.observationElectrodes is None:
            out=[]
            if self.hasInjections():
                if self.hasDipoleInjections() and self.hasDipoleMeasurements():
                    for A,B,M,N in self.tokenIterator():
                        if out.count(M) ==0:
                            out.append(M)
                        if out.count(N) ==0:
                            out.append(N)

                elif not self.hasDipoleInjections() and self.hasDipoleMeasurements():
                    for A,M,N in self.tokenIterator():
                        if out.count(M) ==0:
                            out.append(M)
                        if out.count(N) ==0:
                            out.append(N)
                elif self.hasDipoleInjections() and not self.hasDipoleMeasurements():
                    for A,B, M in self.tokenIterator():
                        if out.count(M) ==0:
                            out.append(M)
                else:
                    for A,M in self.tokenIterator():
                        if out.count(M) ==0:
                            out.append(M)
            else:

                if self.hasDipoleMeasurements():
                    for M,N in self.tokenIterator():
                        if out.count(M) ==0:
                            out.append(M)
                        if out.count(N) ==0:
                            out.append(N)
                else:
                    for M, in self.tokenIterator():
                        if out.count(M) ==0:
                            out.append(M)
            self.observationElectrodes=out
        return self.observationElectrodes   

    def getObservationElectrodeIndex(self, M):
        """
        returns the index of the injection station M in self.getObservationElectrodes()
        """
        l=self.getObservationElectrodes()
        if M in l:
            return l.index(M)
        else:
            raise ValueError("Unknown observation electrode %s"%M)
        
    def getObservations(self, A=None, B=None, insertSource=False):
        """
        returns a list of the observations (M, N) using injection (A,B)
        """
        if isinstance(A, tuple):
            A, B= A[0], A[1]
            
        out=[]
        if self.hasInjections():
            if self.hasDipoleInjections() and self.hasDipoleMeasurements():
                for A1,B1,M,N in self.tokenIterator():
                    if A1 == A and B1 == B:
                        if insertSource:
                            t=(A, B, M, N)
                        else:
                            t=(M,N)
                        if not t in out:
                                out.append(t)

            elif not self.hasDipoleInjections() and self.hasDipoleMeasurements():
                for A1,M,N in self.tokenIterator():
                    if A1 == A:
                        if insertSource:
                            t=(A, M, N)
                        else:
                            t=(M,N)
                        if not t in out:
                            out.append(t)

            elif self.hasDipoleInjections() and not self.hasDipoleMeasurements():
                for A1,B1, M in self.tokenIterator():
                    if A1 == A and B1 == B:
                        if insertSource:
                            t=(A, B, M)
                        else:
                            t=M
                        if not t in out:
                            out.append(t)
            else:
                for A1, M in self.tokenIterator():
                    if A1 == A:
                        if insertSource:
                            t=(A, M)
                        else:
                            t=M
                        if not t in out:
                            out.append(t)
        else:
            if self.hasDipoleMeasurements():
                for M,N in self.tokenIterator():
                        if out.count((M,N)) ==0:
                            out.append((M,N))
            else:
                for M, in self.tokenIterator():
                    if out.count(M) ==0:
                        out.append(M)

        return out   


    def makePrediction(self, values={}, valuesKeysAreStationKeys = True):
        """
        creates prediction for the values[S] predicted by the (single electrode!) injection at station S in getListOfInjectionStations()
        at the measurement stations in order getStationNumeration(), ie. values[S][i] is the predicted potential when injecting at station S and measuring at station getKeyOfStationNumber[i]  
        """
        out={}
        if self.hasInjections():
            if self.hasDipoleInjections() and self.hasDipoleMeasurements():
                for A,B,M,N in self.tokenIterator():
                    if valuesKeysAreStationKeys:
                        iA, iB = A, B
                    else:
                        iA, iB = self.getStationNumber(A), self.getStationNumber(B)
                    out[(A,B,M,N)] = (values[iA][self.getStationNumber(M)]-values[iB][self.getStationNumber(M)]-
                                      values[iA][self.getStationNumber(N)]+values[iB][self.getStationNumber(N)])
            elif not self.hasDipoleInjections() and self.hasDipoleMeasurements():
                for A,M,N in self.tokenIterator():
                    if valuesKeysAreStationKeys:
                        iA =  A
                    else:
                        iA = getStationNumber(A)

                    out[(A,M,N)] = values[iA][self.getStationNumber(M)]-values[iA][self.getStationNumber(N)]
            elif self.hasDipoleInjections() and not self.hasDipoleMeasurements():
                for A,B,M in self.tokenIterator():
                    if valuesKeysAreStationKeys:
                        iA, iB = A, B
                    else:
                        iA, iB = self.getStationNumber(A), self.getStationNumber(B)
                    out[(A,B,M)] = values[iA][self.getStationNumber(M)]-values[iB][self.getStationNumber(M)]
            else:
                for A,M in self.tokenIterator():
                    if valuesKeysAreStationKeys:
                        iA =  A
                    else:
                        iA= getStationNumber(A)
                    out[(A,M)] = values[iA][self.getStationNumber(M)]
        else:
            if self.hasDipoleMeasurements():
                for M,N in self.tokenIterator():
                    out[(M,N)] = values[self.getStationNumber(M)]-values[self.getStationNumber(N)]
            else:
                for M in self.tokenIterator():
                    out[M] = values[self.getStationNumber(M)]
            
        return out
    
        
    def makeResistencePrediction(self, values={}, **kwargs ):
        """
        creates prediction for the resistance based on the potentials values[S] predicted by the (single electrode!) injection at station S in getListOfInjectionStations()
        at the measurement stations in order getStationNumeration(), ie. values[S][i] is the predicted potential when injecting at station S and measuring at station getKeyOfStationNumber[i]  
        """
        return self.makePrediction(values, **kwargs)
    
    def saveDataToCSV(self, fn):
        raise NotImplementedErrior
        f=open(fn, 'w')
        if self.hasDipoleInjections():
            for A, B, M, N in self.getSchedule():
                f.write("%d, %d, %d, %d, %10.15e\n"%(A,B,M,N, data[(A, B, M, N)]))
        else:
            for A, B, M in self.getSchedule():
                f.write("%d, %d, %d, %10.15e\n"%(A,B,M, data[(A, B, M, N)]))
        f.close()  

