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
                           dipoleMeasurements=True, delimiter=',', commend='#', printInfo=True):
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
    data=SurveyData(stations=stations, observations=columns, dipoleInjections=dipoleInjections, dipoleMeasurements=dipoleMeasurements, hasInjections=hasInjections)
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
                F=[ float(ll[i+ns*c]) for i in range(nc) ]
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
            
            T=(A, M, N)

    or if dipoleMeasurements and dipoleInjections are not set 
            
            T=(A, M)
    
    and (columns[0], columns[1], ....) refers to the data entries at each observation.
    
    
    The location of the injecting or observing electrode A is given by stationlocations[A]


    """
    OBSTYPES= ['R', 'E', 'ETA', 'GAMMA', 'ERR_R', 'ERR_E', 'ERR_ETA', 'ERR_GAMMA', 'RELERR_R', 'RELERR_E', 'RELERR_ETA', 'RELERR_GAMMA']
    def __init__(self, stations={}, observations=[], dipoleInjections=True, dipoleMeasurements=True,  hasInjections=True, default_rel_error=0.01):
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
    
    @classmethod 
    def checkObservationType(cls, obs):
        if isinstance(obs, list):
            return all([cls.checkObservationType(o) for o in obs])
        else:
            return obs in cls.OBSTYPES
    
    
    def getStationNumeration(self):
        """
        returns list of the station identifiers in a particular order
        """
        return self.station_key
    
    def getStationNumber(self, k):
        """
        returns the station number from station identfier
        """
        return self.station_index[k]
    
    def getKeyOfStationNumber(self, k):
        return self.station_key[k]
    
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
    def getStationLocation(self, s):
        return self.stations[s]

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
    
    def getResistenceData(self, token):
        return self.getDataRecord(token, datatype='R')

    def getResistenceRelError(self, token):
        if self.hasDataType("RELERR_R"):
            return self.getDataRecord(token, datatype='RELERR_R')
        elif self.hasDataType("ERR_R"):
            r=self.getDataRecord(token, datatype='R')
            e=self.getDataRecord(token, datatype='ERR_R')
            if abs(r)>0:
                return abs(e/r)
            else:    
                return e
        else:
            return self.default_rel_error
            

    def getFieldIntensityData(self, token):
        return self.getDataRecord(token, datatype='E')
    def getFieldData(self, token):
        return self.getDataRecord(token, datatype='E0'), self.getDataRecord(token, datatype='E1'), self.getDataRecord(token, datatype='E2')
    def getGravityData(self, token):
        return self.getDataRecord(token, datatype='gz')

    def getChargeabilityData(self, token):
        return self.getDataRecord(token, datatype='gamma')
    
    def getFieldDifferenceData(self, token):
        return self.getDataRecord(token, datatype='dE0'), self.getDataRecord(token, datatype='dE1'), self.getDataRecord(token, datatype='dE2')
        
    def tokenIterator(self):
        """
        returns a list of tokens (e.g. list of (A,B,M,N)s)
        """
        return self.data.keys()
    
    def injectionIterator(self):
        """
        returns an iterator (in form of list) of the injection dipols or stations:
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
        returns a list of the stations used for injection
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


    def makePrediction(self, values={}):
        """
        creates prediction for the values[S] predicted by the (single electrode!) injection at station S in getListOfInjectionStations()
        at the measurement stations in order getStationNumeration(), ie. values[S][i] is the predicted potential when injecting at station S and measuring at station getKeyOfStationNumber[i]  
        """
        out={}
        if self.hasInjections():
            if self.hasDipoleInjections() and self.hasDipoleMeasurements():
                for A,B,M,N in self.tokenIterator():
                    out[(A,B,M,N)] = values[A][self.getStationNumber(M)]-values[B][self.getStationNumber(M)]-values[A][self.getStationNumber(N)]+values[B][self.getStationNumber(N)]
            elif not self.hasDipoleInjections() and self.hasDipoleMeasurements():
                for A,M,N in self.tokenIterator():
                    out[(A,M,N)] = values[A][self.getStationNumber(M)]-values[A][self.getStationNumber(N)]
            elif self.hasDipoleInjections() and not self.hasDipoleMeasurements():
                for A,B,M in self.tokenIterator():
                    out[(A,B,M)] = values[A][self.getStationNumber(M)]-values[B][self.getStationNumber(M)]
            else:
                for A,M in self.tokenIterator():
                    out[(A,M)] = values[A][self.getStationNumber(M)]
        else:
            if self.hasDipoleMeasurements():
                for M,N in self.tokenIterator():
                    out[(M,N)] = values[self.getStationNumber(M)]-values[self.getStationNumber(N)]
            else:
                for M in self.tokenIterator():
                    out[M] = values[self.getStationNumber(M)]
            
        return out
    
        
    def makeResistencePrediction(self, values={}):
        """
        creates prediction for the resistance based on the potentials values[S] predicted by the (single electrode!) injection at station S in getListOfInjectionStations()
        at the measurement stations in order getStationNumeration(), ie. values[S][i] is the predicted potential when injecting at station S and measuring at station getKeyOfStationNumber[i]  
        """
        return self.makePrediction(values)
    
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

