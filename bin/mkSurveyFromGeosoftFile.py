#!/usr/bin/python3
"""
Fingal - creates configuration, station and data file from Geosoft file. 
This is for the field intensity data only. 

by l.gross@uq.edu.au, Dec 2020.
"""
import argparse
import numpy as np
from pathlib import Path
import sys,os 
from fingal import FindNearestElectrode

sys.path.append(os.getcwd())

parser = argparse.ArgumentParser(description='creates survey data sets (station file and data file) for a FullWaver survey from geosoft survey data', epilog="version 1/8/2021")
parser.add_argument('--resmin', '-r', dest='resmin', metavar='RESmin', type=float, default=0.0, help='ignore records with resitivity lower then Resmin')
parser.add_argument('--fullwaver', '-f', dest='fullwaver', action='store_true', default=False, help='create fullwaver data set')
parser.add_argument('--xyz', '-x', dest='directional',  action='store_true', default=False, help='create directional quantities for Fullwaver')
parser.add_argument('--nocenter', '-n', dest='nocenter',  action='store_true', default=False, help='electrode locations are centered to zero.')
parser.add_argument('--datafile', '-d',dest='datafile', metavar='DATAFILE', type=str, help='data file name to read from. if not given [PROJECT].dat is used.')
parser.add_argument(dest='project', metavar='PROJECT', type=str, help='project name')
#parser.add_argument('--novtx', '-nv', dest='ignorevtx', action='store_true', help="if set vtx points are ignored.")


args = parser.parse_args() 

# colums indices:
AX, AY, AZ =0, 1, 2
BX, BY, BZ =3, 4, 5
MX, MY, MZ =6, 7, 8
NX, NY, NZ =9, 10, 11
I, VP, RES, ETA =26, 27, 28, 32

if args.datafile:
    datafile=args.datafile
else:
    datafile=args.project+".dat"
  
electrodesFN=args.project+".loc"
dataFN=args.project+".csv"
configFN=args.project+".py"
meshfile=args.project+".fly"



    
print("data file %s is opened."%datafile)
TABS=None
getTABS=False
l=[]
datatable=[]
A, B, M, N, Ic, Rho, Eta =[], [], [],  [], [], [] , []

with open(datafile, 'r') as f:
    n=1
    line=f.readline()
    while line:
        if getTABS and TABS is None:
            TABS=l2.split()
            print("Columns:Electrode A Coordinates: ",TABS[AX], TABS[AY], TABS[AZ])
            print("\tElectrode A Coordinates: ",TABS[BX], TABS[BY], TABS[BZ])
            print("\tElectrode B Coordinates: ",TABS[MX], TABS[MY], TABS[MZ])
            print("\tElectrode M Coordinates: ",TABS[NX], TABS[NY], TABS[NZ])
            print("\tInjected current:  ",TABS[I])
            print("\tVP:  ",TABS[VP])
            print("\tResistivity:  ",TABS[RES])
            print("\tChargeability:  ",TABS[ETA])

        l, l2=line.strip(), l
        if l[0].isdigit() or l[0] in ["+", "-", "."]:
            getTABS=True
            d=[ float(f) for f in l.split() ]
            if d[RES] > args.resmin :
                A.append([d[AX], d[AY], d[AZ]])
                B.append([d[BX], d[BY], d[BZ]])
                M.append([d[MX], d[MY], d[MZ]])
                N.append([d[NX], d[NY], d[NZ]])
                Ic.append(d[I])
                Rho.append(d[RES])
                Eta.append(d[ETA])
            else:
                print("line %s ignored as not rho = %s > %s=RESmin."%(n, d[RES], args.resmin))

        n+=1
        line=f.readline()
print("%s records found"%len(A))
A=np.array(A)
B=np.array(B)
M=np.array(M)
N=np.array(N)
Ic=np.array(Ic)
Rho=np.array(Rho)
Eta=np.array(Eta)
print("resistivity range = %e - %e"%(Rho.min(), Rho.max()))
rho_bar=1./((1./Rho).mean())
print("mean resistivity = %e"%(rho_bar))
L=np.sqrt(A[:,0]**2+A[:,1]**2+A[:,2]**2).max()
print("characteristic length = ",L)



electrodes={} # dictionary of found electrodes: id -> coordinates
surveys={}    # surveys[(eidA, eidB)][(eidM, eidN)] = (sum Rho, sum Eta, counter) 


values={}
ne=0
ff=[]
nn=0

for s in range(len(A)):
    if abs(Ic[s]) >0:
        eidA, edist=FindNearestElectrode(A[s,0], A[s,1], A[s,2], electrodes)
        if edist > 1.e-8*L:
            ne+=1
            electrodes[ne]=(A[s,0], A[s,1], A[s,2])
            eidA=ne

        eidB, edist=FindNearestElectrode(B[s,0], B[s,1], B[s,2], electrodes)
        if edist > 1.e-8*L:
            ne+=1
            electrodes[ne]=(B[s,0], B[s,1], B[s,2])
            eidB=ne

        eidM, edist=FindNearestElectrode(M[s,0], M[s,1], M[s,2], electrodes)
        if edist > 1.e-8*L:
            ne+=1
            electrodes[ne]=(M[s,0], M[s,1], M[s,2])
            eidM=ne

        eidN, edist=FindNearestElectrode(N[s,0], N[s,1], N[s,2], electrodes)
        if edist > 1.e-8*L:
            ne+=1
            electrodes[ne]=(N[s,0], N[s,1], N[s,2])
            eidN=ne
        values[(eidA, eidB, eidM, eidN)] = Rho[s]
        #print (eidA, eidB, eidM, eidN), M[s,0], M[s,1], M[s,2], ":", N[s,0], N[s,1], N[s,2]

        if not (eidA, eidB) in surveys.keys():
            surveys[(eidA, eidB)]={}
        v=Rho[s]
        if (eidM, eidN) in surveys[(eidA, eidB)].keys():
            surveys[(eidA, eidB)][(eidM, eidN)] = ( surveys[(eidA, eidB)][(eidM, eidN)][0]+Rho[s], surveys[(eidA, eidB)][(eidM, eidN)][1]+Eta[s], surveys[(eidA, eidB)][(eidM, eidN)][2]+1)
        else:
            surveys[(eidA, eidB)][(eidM, eidN)] = (Rho[s], Eta[s], 1)
        nn+=1

# apply averaging at each data record:
for inject in surveys:
    for obs in surveys[inject]:
        rho1=surveys[inject][obs][0]/surveys[inject][obs][2]
        m1=surveys[inject][obs][1]/surveys[inject][obs][2]/1000.
        surveys[inject][obs] = (rho1, m1)

if not args.nocenter :
    n=len(electrodes)
    x0=sum( [ electrodes[ie][0] for ie in electrodes  ])/n
    y0=sum( [ electrodes[ie][1] for ie in electrodes  ])/n
    print("center at (%s, %s)"%(x0, y0))
    
    for ie in electrodes:
        electrodes[ie]=(electrodes[ie][0]-x0, electrodes[ie][1]-y0, electrodes[ie][2])
print("%s records found "%nn)
with open(electrodesFN, 'w') as f:
    for ie in electrodes:
        f.write("%d, %g, %g, %g\n"%(ie, electrodes[ie][0], electrodes[ie][1], electrodes[ie][2]))
print("%s electrode locations written to file %s."%(len(electrodes), electrodesFN) )

if args.fullwaver and not args.directional:
                    
    absEs={}
    absGs={}
    for surv in surveys:
        absEs[surv]={}
        absGs[surv]={}
        ee=[]
        ee2=[]
        for m,n in surveys[surv]:
            if ee.count(m)==0:
                ee.append(m)
            else:
                ee2.append(m)
            if ee.count(n)==0:
                ee.append(n)
            else:
                ee2.append(n)
        ee2.sort()
        
        for e in ee2:
            p1=None
            p2=None
            for m,n in surveys[surv]:
                if m==e:
                    p1=n
                if n==e:
                    p2=m
            xe=np.array(electrodes[e])
            xa=np.array(electrodes[surv[0]])
            xb=np.array(electrodes[surv[1]])
            rae=np.linalg.norm(xa-xe)
            rbe=np.linalg.norm(xb-xe)
        
            rho1=surveys[surv][(e,p1)][0]
            m1=surveys[surv][(e,p1)][1]
            x1=np.array(electrodes[p1])
            l1 =np.linalg.norm(xe-x1)
            ra1=np.linalg.norm(xa-x1)
            rb1=np.linalg.norm(xb-x1)
            d1=1./(1./ra1 -1./rb1 -1./rae +1./rbe)
            E1=-rho1/(l1*d1)/2./np.pi
        
            rho2=surveys[surv][(p2,e)][0]
            m2=surveys[surv][(p2,e)][1]
            x2=np.array(electrodes[p2])
            l2 =np.linalg.norm(xe-x2)
            ra2=np.linalg.norm(xa-x2)
            rb2=np.linalg.norm(xb-x2)
            d2=1./(1./ra2 -1./rb2 -1./rae +1./rbe)
            E2=-rho2/(l2*d2)/2./np.pi
 
            gamma1=m1/(1-m1)
            gamma2=m2/(1-m2)
            E=(E1**2+E2**2)**0.5
            gamma=E1**2/(E1**2+E2**2)*gamma1+E2**2/(E1**2+E2**2)*gamma2
            if abs(E1)>0 and abs(E2)>0:
                #print e,"->", p1, rho1, m1, l1, p2, rho2, m2, l2," -> ",E, M
                absEs[surv][e]=E
                absGs[surv][e]=gamma

    FMT=("%d, "*3+"%g, %g\n")
    with open(dataFN,'w') as f:
        nout=0
        for a,b in absEs:
            for m in absEs[(a,b)]:
                f.write(FMT%(a, b, m, absEs[(a,b)][m],absGs[(a,b)][m] ))
                #f.write(FMT%(electrodes[a][0],electrodes[a][1],electrodes[a][2], electrodes[b][0],electrodes[b][1],electrodes[b][2], electrodes[m][0],electrodes[m][1],electrodes[m][2], absEs[(a,b)][m],absGs[(a,b)][m] ))
                nout+=1
        print("%d data records written to %s."%(nout, dataFN))
    
else:
    print("output not supported!")

with open(configFN, 'w') as fout:
    fout.write("meshfile = '%s'\n"%meshfile)
    fout.write("stationfile = '%s'\n"%electrodesFN)
    fout.write("stationdelimiter = ','\n")
    fout.write("usesStationCoordinates = False\n")
    fout.write("stationsFMT = None\n")
    fout.write("dipoleInjections = True\n")
    fout.write("dipoleMeasurements = False\n")
    fout.write("datadelimiter = ','\n")
    fout.write("datafile = '%s'\n"%dataFN)
    fout.write("schedulefile = '%s'\n"%dataFN)
    fout.write("restartfile = 'restarter'\n")


    fout.write("sigma0=%s\n"%(1/rho_bar))
    fout.write("sigma_background=sigma0\n")
    fout.write("eta0=0.01\n")
    fout.write("eta_background=0.\n")


    fout.write("L_stations=1.\n")
    fout.write("adjustStationLocationsToElementCenter=True\n")
    fout.write("useLogDefect=True\n")
    fout.write("dipoles=False\n")
    fout.write("gamma0=eta0/(1-eta0)\n")
    fout.write("gamma_background=eta_background/(1-eta_background)\n")
    fout.write("tolerance=1e-4\n")
    fout.write("interpolation_order=3\n")
    fout.write("w1=0\n")
    fout.write("alpha0=1.\n")
    fout.write("alpha1=(130.)**2\n")
    fout.write("alpha0gamma=alpha0\n")
    fout.write("alpha1gamma=alpha1\n")
    fout.write("w1gamma=w1\n")

    fout.write("w0=1e-10\n")    
    fout.write("w0gamma=w0\n")

    fout.write('outfile="sigma"\n')
