#!/usr/bin/python3
from esys.pycad import *
from esys.pycad.gmsh import Design
import os
from math import *
from esys.finley import MakeDomain, ReadMesh, ReadGmsh 
from esys.escript import gmshGeo2Msh, length, Function
from esys.weipa import saveSilo
import numpy as np
import argparse
desc="""
extracting the surface elements of the CORE mesh given in fly format and
generates a gmsh GEO file to mesh an region outside of the core. By default
the GEO sets up a bounding box around the CORE but the file can be edited to
add more features (eg. layering). When the outer region has been meshed the
new mesh can be merged with CORE mesh.
"""
parser = argparse.ArgumentParser(description=desc, epilog="l.gross@uq.edu.au version 21/12/2020")
parser.add_argument(dest='coremeshfile', metavar='CORE', type=str, help='core region in fly format (input)')
parser.add_argument(dest='geofile', metavar='GEO', type=str, help='geofile file be generated. (output)')
parser.add_argument('--sanitycheck', '-s',  dest='sanitycheck', action='store_true', default=False, help="test for voids.")
parser.add_argument('--tolerance', '-t',  dest='tolerance',  default=1e-8, help="test tolerance")

args = parser.parse_args()

ff=open(args.coremeshfile,'r')
ff.readline()
# read nodes:
numNodes=int(ff.readline().strip().split(" ")[1])
nodes={}
for n in xrange(numNodes):
    node_id, a, tag, x,y,z = ff.readline().strip().split(" ")
    nodes[int(node_id)]=np.array([float(x), float(y), float(z)])

# get the bounding box:    
X_min=min( [ C[0] for C in nodes.itervalues() ] )
X_max=max( [ C[0] for C in nodes.itervalues() ] )

Y_min=min( [ C[1] for C in nodes.itervalues() ] )
Y_max=max( [ C[1] for C in nodes.itervalues() ] )

Z_min=min( [ C[2] for C in nodes.itervalues() ] )
Z_max=max( [ C[2] for C in nodes.itervalues() ] )
DIAM=max(X_max-X_min, Y_max-Y_min, Z_max-Z_min)



print("file "+args.coremeshfile+" opened.")
print("number of input nodes = ",numNodes)
print("Bounding box of input domain:")
print("X range =",X_min, X_max)
print("Y range =",Y_min, Y_max)
print("Z range =",Z_min, Z_max)
print("Domain diameter =",DIAM)


TETFACES=[ (0,2,1) , (1,2,3), (0,1,3), (3,2,0)]
tagmax=0
# jump over the interior elements:
elementfaces={}
numElements=int(ff.readline().strip().split(" ")[1])
for n in xrange(numElements):
    eid, tag, n1,n2,n3, n4 = ff.readline().strip().split(" ")
    nn=[int(n1), int(n2), int(n3), int(n4)]
    for f in TETFACES:
        m=[ nn[f[0]], nn[f[1]], nn[f[2]] ]
        m.sort()
        m=tuple(m)
        if elementfaces.has_key(m):
            elementfaces[m]=None
        else:
            elementfaces[m]=tuple([ nn[f[0]], nn[f[1]], nn[f[2]] ])
    tagmax=max(tagmax,int(tag) )
print(numElements," interior elements found.")

facepatches=[ elementfaces[e] for e in elementfaces if elementfaces[e] is not None ]
print(len(facepatches)," surface patches found.")
del elementfaces
# read surface elements:
numSurfaceElements=int(ff.readline().strip().split(" ")[1])
for n in xrange(numSurfaceElements):
    eid, tag, n1,n2,n3 = ff.readline().strip().split(" ")
    tagmax=max(tagmax,int(tag) )
print("maximal tag used = ", tagmax)

if args.sanitycheck:
    # check for holes:
    FAILED=False
    for fp in facepatches:
        # check area:
        X=nodes[fp[1]]-nodes[fp[0]]
        Y=nodes[fp[2]]-nodes[fp[0]]
        
        A2=(X[1]*Y[2]-X[2]*Y[1])**2+(X[2]*Y[0]-X[0]*Y[2])**2+(X[0]*Y[1]-X[1]*Y[0])**2
        LX=np.linalg.norm(X)
        LY=np.linalg.norm(Y)
        if min(LX,LY) < args.tolerance*DIAM:
            print("zero element edge detected near ",(nodes[fp[1]]+nodes[fp[0]]+nodes[fp[2]])/3)
            FAILED=True
        if A2 < args.tolerance*LX*LX*LY*LY:
            print("zero a area trangle detected near ",(nodes[fp[1]]+nodes[fp[0]]+nodes[fp[2]])/3)
            FAILED=True
        # add to conenctivity
    if FAILED:
        print("WARNING: triangle edge length and area test failed.")
    else:
        print("triangle edge length and area test passed.")
    # === check for holes by node connectvity :
    comps=[]
    for fp in facepatches:
        fps=set(fp)
        c0=None
        for c in comps:
            if not c.isdisjoint(fps):
                c0=c
                break
        if c0 is None:
            comps.append(fps)
        else:
            c0|=fps

    sp=0
    while sp+1 < len(comps):
        for p in range(sp+1, len(comps) ):
            if not comps[sp].isdisjoint(comps[p]):
                comps[sp]|=comps[p]
                sp=0
                del comps[p]
                break
    if len(comps)<2:
        print("single surface component detected! All good")
    else:
        print("WARNING: ",len(comps)," surface components are found.")
        comps.sort(key=lambda o : len(o))
        for i in range(len(comps)):
            X=np.zeros((3,))
            N=0
            for p in comps[i]:
                    X+=nodes[p]
                    N+=1
            print(">> component ",i," with ",len(comps[i])," edges. center near ",X/N)
    del comps
    # ========= check for holes by edge connectvity :
    comps=[]
    for fp in facepatches:
        fps=set( [ ( min(fp[i1], fp[(i1+1)%3]),  max(fp[i1], fp[(i1+1)%3]) ) for i1 in xrange(3) ] )
        c0=None
        for c in comps:
            if not c.isdisjoint(fps):
                c0=c
                break
        if c0 is None:
            comps.append(fps)
        else:
            c0|=fps

    sp=0
    while sp+1 < len(comps):
        for p in range(sp+1, len(comps) ):
            if not comps[sp].isdisjoint(comps[p]):
                comps[sp]|=comps[p]
                sp=-1
                del comps[p]
                break
        sp+=1
    if len(comps)<2:
        print("single surface component based on edge connectivity detected ! All good")
    else:
        print("WARNING: ",len(comps)," surface  based on edge connectivity components are found.")
        comps.sort(key=lambda o : len(o))
        for i in range(len(comps)):
            X=np.zeros((3,))
            N=0
            for p in comps[i]:
                    X+=nodes[p[0]]+nodes[p[1]]
                    N+=2
            print(">> component ",i," with ",len(comps[i])," nodes. center at ",X/N)
    del comps
    
    # now we check if all edges are contained in exactly two triangles:
    FAILED=False
    lines={}
    for e in facepatches:
        for i1 in xrange(3):
            i2=(i1+1)%3
            e0=min(e[i1], e[i2])
            e1=max(e[i1], e[i2])
            if not lines.has_key((e0, e1)):
                lines[(e0, e1)]=1
            else:
                lines[(e0, e1)]+=1
    for e in lines:
        if not lines[e] == 2:
            FAILED=True
            print("WARNING: edge ",e, " is shared by ",lines[e]," face triangles near ",(nodes[e[1]]+nodes[e[0]])/2)
    
    if FAILED:
        print("FAILED: edge counter test failed.")
    else:
        print("edge counter test passed.")

    #print("starting checking intercept of edges. Please be patient.")
    #FAILED=False
    #DIST=DIAM
    #WHERE=None
    
    
    #for e in lines:
        #X=nodes[e[0]]
        #V=nodes[e[1]]-nodes[e[0]]
        #ATOL=(V[0]**2+V[0]**2+V[0]**2)*args.tolerance

        #for f in lines:
            #if f[0]<e[0] and not ( f[0]==e[0] or f[0]==e[1] or f[1]==e[0] or f[1]==e[1] ):
                    #Y=nodes[f[0]]
                    #W=nodes[f[1]]-nodes[f[0]]

                    #A=np.vstack([V,-W])
                    #try :
                        #t=np.linalg.solve(np.dot(A,A.T), np.dot(A,Y-X) )
                        #if 0<=t[0] and t[0]<=1. and 0<=t[1] and t[1]<=1. :
                            #d=numpy.linalg.norm(np.dot(A.T,t)-(Y-X))
                            #if d < DIST:
                                #DIST=d
                                #WHERE=X+t[0]*V
                                #print("edge distance reduced to ",DIST," at ", WHERE, t)
                            #if d*d < ATOL:
                                #print "intercept at ",X+t[0]*V, "debug d, tX, tY =",d, t
                                #FAILED=True
                    #except np.linalg.linalg.LinAlgError:
                        #pass

                    ##t, d, r, s = np.linalg.lstsq(np.vstack([V,-W]).T, Y-X)
                    ###print t, d, r, s 
                    ##if r == 2:
                        ##if 0<=t[0] and t[0]<=1. and 0<=t[1] and t[1]<=1. :
                            ##d=d[0]
                            ##if d < DIST:
                                ##DIST=d
                                ##WHERE=X+t[0]*V
                                ##print "edge distance reduced to ",DIST," at ", WHERE, t, s
                            ##if d < ATOL:
                                ##print "intercept at ",X+tX*V, "debug d, tX, tY =",d, t
                                ##FAILED=True
                    ##elif r==1:
                        ##print t, d, r, s
                    ##else:
                        ##print "length zero edge at ",X
                        ##FAILED=True
                        
                    ##D=-V[0]*W[1]+V[1]*W[0]

                    ##if abs(D)>ATOL: 
                        ##tX=((X[0]-Y[0])*W[1]-W[0]*(X[1]-Y[1]))/D
                        ##if 0<=tX and tX<=1.:
                            ##tY=(-V[0]*(X[1]-Y[1])+V[1]*(X[0]-Y[0]))/D
                            ##if 0<=tY and tY<=1.:
                                ##d=abs(X[2]-Y[2]+tX*V[2]-tY*W[2])
                                ###print f,e, tX, tY, d
                                ##if d < DIST:
                                    ##DIST=d
                                    ##WHERE=X+tX*V
                                    ##print "edge distance reduced to ",DIST," at ", WHERE
                                ##if d*d < ATOL:
                                    ##print "intercept at ",X+tX*V, "debug d, tX, tY =",d, tX, tY
                                    ##FAILED=True
    #del lines
    #if WHERE is not None:
        #print "distance closest segment = ",DIST," at ",WHERE
    #if FAILED:
        #print "FAILED: edge intercept test failed."
    #else:
        #print " edge intercept test passed."


points={}
lines={}
text=""

ns=1
np=1
nl=1
for e in facepatches:
    # register points
    for i in xrange(3):
        if not points.has_key(e[i]):
            points[e[i]]=np
            np+=1
    l=[0,0,0]
    s=["", "" , ""]
    # register lines
    for i1 in xrange(3):
        i2=(i1+1)%3
        e0=min(e[i1], e[i2])
        e1=max(e[i1], e[i2])
        if not lines.has_key((e0, e1)):
            lines[(e0, e1)]=nl
            nl+=1
        l[i1]=lines[(e0, e1)]  
        if e[i1] > e[i2]:
            s[i1]="-"

    text+="ll=newreg ;Line Loop(ll) = {%sl%s, %sl%s, %sl%s};\n"%(s[0], l[0],  s[1], l[1], s[2], l[2])
    text+="coresurfaces[%s] = newreg; Plane Surface(coresurfaces[%s]) = {-ll};\n"%(ns-1,ns-1)
    ns+=1
# new we create the points:
text0="// start of defining the surface loops of the core \n"
np=0
for p in points:
    text0+="p%s = newp; Point(p%s) = {%10.15e,  %10.15e,  %10.15e};\n"%(points[p], points[p], nodes[p][0], nodes[p][1], nodes[p][2])  
    np+=1
print(np," points generated.")
nl=0
for l in lines:
    text0+="l%s = newl; Line(l%s) = {p%s, p%s}; Transfinite Line {l%s} = 2;\n"%(lines[l], lines[l], points[l[0]], points[l[1]],  lines[l])      
    nl+=1
    
print(nl-1," lines generated.")
print(ns-1," surface loops generated.")
    
# now we are ready to add the pedding region:
text0+=text
text0+="coresl=newreg; Surface Loop(coresl) = {coresurfaces[]};\n"

out ="CoreXmin=%e;\n"%X_min
out+="CoreXmax=%e;\n"%X_max
out+="CoreLX=CoreXmax-CoreXmin;\n"
out+="CoreYmin=%e;\n"%Y_min
out+="CoreYmax=%e;\n"%Y_max
out+="CoreLY=CoreYmax-CoreYmin;\n"
out+="CoreZmin=%e;\n"%Z_min
out+="CoreZmax=%e;\n"%Z_max
out+="CoreLZ=CoreZmax-CoreZmin;\n"
out+="CoreDiag=%e;\n"%DIAM
out+="""
// This is an example of a padding box around the core.  please modify
BBXmin=CoreXmin-CoreLX/2;
BBXmax=CoreXmax+CoreLX/2;
BBYmin=CoreYmin-CoreLY/2;
BBYmax=CoreYmax+CoreLY/2;
BBZmin=CoreZmin-CoreLZ/2;
BBZmax=CoreZmax+CoreLZ/2;
BBmeshSize=CoreLX/5;

Point(1) = {BBXmin, BBYmin, BBZmin, BBmeshSize};
Point(2) = {BBXmax, BBYmin, BBZmin, BBmeshSize};
Point(3) = {BBXmin, BBYmax, BBZmin, BBmeshSize};
Point(4) = {BBXmax, BBYmax, BBZmin, BBmeshSize};
Point(5) = {BBXmin, BBYmin, BBZmax, BBmeshSize};
Point(6) = {BBXmax, BBYmin, BBZmax, BBmeshSize};
Point(7) = {BBXmin, BBYmax, BBZmax, BBmeshSize};
Point(8) = {BBXmax, BBYmax, BBZmax, BBmeshSize};
Line(9) = {1, 2};
Line(10) = {2, 4};
Line(11) = {4, 3};
Line(12) = {3, 1};
Line(13) = {5, 6};
Line(14) = {6, 8};
Line(15) = {8, 7};
Line(16) = {7, 5};
Line(17) = {1, 5};
Line(18) = {2, 6};
Line(19) = {3, 7};
Line(20) = {4, 8};

Line Loop(21) = {9, 10, 11, 12};
Plane Surface(22) = {-21};
Physical Surface(\"bottom\") = {22};

Line Loop(23) = {13, 14, 15, 16};
Plane Surface(24) = {23};
Physical Surface(\"top\") = {24};

Line Loop(25) = {9, 18, -13, -17};
Plane Surface(26) = {25};
Physical Surface(\"front\") = {26};

Line Loop(27) = {-11, 20, 15, -19};
Plane Surface(28) = {-27};
Physical Surface(\"back\") = {28};

Line Loop(29) = {17, -16, -19, 12};
Plane Surface(30) = {29};
Physical Surface(\"left\") = {30};

Line Loop(31) = {18, 14, -20, -10};
Plane Surface(32) = {-31};
Physical Surface(\"right\") = {32};

// surface loop of the outer boundary of the domain  
Surface Loop(33) = {22, 24, 26, 28, 30, 32};

"""
out+=text0
out+="""
// generate volume without core:
paddingVol=newv ;
Volume(paddingVol) = {33, coresl};
Physical Volume(\"padding\") = {paddingVol};
"""

open(args.geofile,'w').write(out)
print("core surface geometry written to "+args.geofile)
