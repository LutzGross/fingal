#!/usr/bin/python3
#
#  This is a simple COMSOL MPH to escript fly file converter
#
import os, sys
from fingal import makeTagMap
from esys.finley import ReadMesh 
from esys.weipa import saveSilo
from esys.escript import ReducedFunction
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='COMSOL MPH text file (Mesh object) to esys-escript fly file. Currently Tets and Tris are supported only', epilog="l.gross@uq.edu.au 11/12/20")
parser.add_argument('mph', metavar='mphtxt_file', type=str,
                    help='name of MPH text input file')
parser.add_argument('fly', metavar='fly_file', type=str,
                    help='name FLY out file')
parser.add_argument('--pointsfile', '-pf', dest='points', metavar='pfile', type=str, default=None,
                    help='if set the name of file to save location and tags of point sources with row format s, x, y, z')
parser.add_argument('--novtx', '-nv', dest='ignorevtx', action='store_true', help="if set vtx points are ignored when creating the fly file.")
parser.add_argument('--silo', '-s',  dest='silofile', metavar='SILO', help="silo file for saving mesh file for visualization (no extension) (output if set).")

args = parser.parse_args()
fly=open(args.fly, "w")
mph=open(args.mph, "r")
print("open COMSOL file ",args.mph)
######
RECORD_MESH="4 Mesh"

#  jump over the header:
line=mph.readline()
if line.startswith("#"):  # save the header
    fly.write(line)
else:
    fly.write("Created by COMSOL \n")

while not line.startswith(RECORD_MESH):
    line=mph.readline()
print("mesh object found.")
version=int(mph.readline().split("#")[0])
sdim=int(mph.readline().split("#")[0])
numNodes=int(mph.readline().split("#")[0])
minNodeID=int(mph.readline().split("#")[0])

# find the next line that that is not empty and does not start with a "#":
line=mph.readline().strip()
while len(line)==0 or line.startswith("#"):
        line=mph.readline().strip()
# collect all coordinates in a string
nodestr=""
n=0
while len(line)>0:
    nodestr+=line.split("#")[0]+ " "
    n+=1
    line=mph.readline().strip()

nodesX=np.fromstring(nodestr, dtype=float, sep=' ')
assert numNodes*sdim == len(nodesX), "number of nodes %s and spatial dimension %s specified in header do not match number of node coordinates values read %s."%(numNodes, sdim, len(nodes))
nodesX=nodesX.reshape((numNodes, sdim))
# write to fly:

fly.write("%1.1dD-Nodes %d\n"%(sdim, numNodes))
for k in range(numNodes):
    if sdim == 2:
        fly.write("%d %d %d %10.15e %10.15e\n"%(k+minNodeID, k+minNodeID, 0, nodesX[k,0], nodesX[k,1]))
    else:
        fly.write("%d %d %d %10.15e %10.15e %10.15e\n"%(k+minNodeID, k+minNodeID, 0, nodesX[k,0], nodesX[k,1], nodesX[k,2]))
print(numNodes, " nodes were written.")
# types 
numTypes=int(mph.readline().split("#")[0])

print("number of element types found is ",numTypes)

diracpoints={}
interiorelements={}
faceelements={}
contactelements={}
contactelements['n']=0

typecount=0
while typecount < numTypes:
    # find the first line that is not empty and is not starting with #  
    mph.readline().strip()
    while len(line)==0 or line.startswith("#"):
        line=mph.readline().strip()
    elmtype=line.split("#")[0].split(" ")[1]

    line=mph.readline().strip()
    while len(line)==0 or line.startswith("#"):
        line=mph.readline().strip()
    numNodes=int(line.split("#")[0])
    numElem=int(mph.readline().split("#")[0])
    print(elmtype, "found.")
    
    # collect all element in a string
    nodestr=""
    line=mph.readline().strip()
    while len(line)>0:
        nodestr+=line.split("#")[0]+ " "
        line=mph.readline().strip()

    nodes=np.fromstring(nodestr, dtype=int, sep=' ')
    assert numNodes*numElem == len(nodes), "number of nodes %s and of elements %s specified in header do not match number of node ids read %s."%(numNodes, numElem, len(nodes))
    nodes=nodes.reshape((numElem,numNodes))
    print(numElem, " elements found.")
    
    numTags=int(mph.readline().split("#")[0])
    # collect all tags in a string
    tagstr=""
    line=mph.readline().strip()
    while len(line)>0:
        tagstr+=line.split("#")[0]+ " "
        line=mph.readline().strip()
    tags=np.fromstring(tagstr, dtype=int, sep=' ')
    assert numTags == numElem, "number of tags does not match number of elements."
    assert numTags == len(tags), "number of tags found does not match number of specified tags."%(numNodes, numElem, len(nodes))       

    if elmtype == "vtx":
        if len(diracpoints)>0:
            raise KeyError("Dirac elements already defined.")
        diracpoints['nodes']=nodes
        diracpoints['tags']=tags
        diracpoints['type']='Point1'
        diracpoints['n']=numElem
        
    elif elmtype == "edg":
        if sdim==2:
            if len(faceelements)>0:
                raise KeyError("face elements already defined.")
            faceelements['nodes']=nodes
            faceelements['tags']=tags
            faceelements['type']='Line2'
            faceelements['n']=numElem
            contactelements['type']="Line2_Contact"
    elif elmtype == "tri":
        if sdim==2:
            if len(interiorelements)>0:
                raise KeyError("face elements already defined.")
            interiorelements['nodes']=nodes
            interiorelements['tags']=tags
            interiorelements['type']='Tri3'
            interiorelements['n']=numElem
        else:
            if len(faceelements)>0:
                raise KeyError("interior elements already defined.")
            faceelements['nodes']=nodes
            faceelements['tags']=tags
            faceelements['type']='Tri3'
            faceelements['n']=numElem
            contactelements['type']="Tri3_Contact"
    elif elmtype == "tet":
        if sdim==3:
            if len(interiorelements)>0:
                raise KeyError("interior elements already defined.")
            interiorelements['nodes']=nodes
            interiorelements['tags']=tags
            interiorelements['type']='Tet4'
            interiorelements['n']=numElem
        else:
            raise KeyError("not sure what to do with tet elements in a 2D mesh.")
    else:
         raise KeyError("Unknown element type "+elmtype)
    typecount+=1

# write it to fly file:
elid=0
for k,v in enumerate([interiorelements, faceelements, contactelements, diracpoints]):
    if k<3 or not args.ignorevtx:
        fly.write("%s %d\n"%(v['type'], v['n']))
        for j in range(v['n']):
            fly.write("%d %d"%(elid+j, v['tags'][j]))
            for n in v['nodes'][j]:
                fly.write(" %d"%n)
            fly.write("\n")
    else:
        print("vtx element ingnored.")
        fly.write("%s %d\n"%(v['type'], 0))
    elid+=v['n']
fly.write("Tags\n")
fly.close()
mph.close()
print("escript file %s generated."%(args.fly))

if args.points:
    nodes=diracpoints['nodes']
    tags=diracpoints['tags']

    f=open(args.points, 'w')
    for k in range(len(tags)):
        idx=nodes[k][0]
        x=nodesX[idx,:]
        if sdim == 2:
            f.write("%d %g %g %g\n"%(tags[k],  x[0], x[1], 0.))
        else:
            f.write("%d %g %g %g\n"%(tags[k],  x[0], x[1], x[2]))
            
    f.close()
    print("source points have been written to "+args.points)


if args.silofile is not None:
    mesh=ReadMesh(args.fly)
    saveSilo(args.silofile,tags=makeTagMap(ReducedFunction(mesh)))
    print(args.silofile+".silo with tags has been generated.")

if args.points:
        print("***WARNING: in COMSOL tags/ids for geometrical node IDs for sources are shown with an offset by 1 to the qIDs in the file. ")
