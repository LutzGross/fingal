#!/usr/bin/python3
import os
import json

from rich import constrain
from sympy.unify.usympy import construct

LOGDIR = "logs/"
FLYDIR = "meshes/"

loglist=[]
for logf in os.listdir(LOGDIR):
    with open(os.path.join(LOGDIR,logf), "r") as f:
        #print(f.readlines()[-1])
        #print(logf, f.read().splitlines()[-1])
        ll = f.readlines()
        loglist.append(json.loads(ll[-1]))
print(len(loglist), " logs read.")
print(loglist)
# get list f files:
meshes =  list( set( [ rec['mesh'] for rec in loglist ] ))
meshes.sort(key=lambda x: int(x[10:x.index(".")]))
contrasts = list(set([ rec['contrast'] for rec in loglist ] ))
contrasts.sort()
ratios = list(set( [ rec['ratio'] for rec in loglist ] ))
ratios.sort()
meshes_NN = {}
for f in meshes:
    with open(f, 'r') as fp:
        lines = fp.readlines()
        for row in lines:
            if row.startswith('3D-Nodes'):
                meshes_NN[f]=int(row[9:])
                break
for c in contrasts[1:]:
    out = "\\toprule\n"
    out += " # nodes"
    for r in ratios:
        out += " & " + str(r)
    out+= "\\\\\n"
    out += "\\midrule\n"
    for f in meshes:
        out+= f"{meshes_NN[f]} "
        for r in ratios:
            out += f" & "
            for rec in loglist:
                if rec['mesh'] == f and rec['contrast'] == c and rec['ratio'] == r:
                    out+= f"{rec['iterations']}"
                    break
        out += " \\\\\n"
    out+="\\bottomrule\n"
    print(out)

for c in contrasts[1:]:
    out = "\\toprule\n"
    out += " # nodes"
    for r in ratios:
        out += " & " + str(r)
    out+= "\\\\\n"
    out += "\\midrule\n"
    for f in meshes:
        out += f"{meshes_NN[f]} "
        for r in ratios:
            out += f" & "
            for rec in loglist:
                if rec['mesh'] == f and rec['contrast'] == c and rec['ratio'] == r:
                    out+= f"{rec['timing_total']:.3g}"
                    break
        out += " \\\\\n"
    out+="\\bottomrule\n"
    print(out)



print(contrasts)
print(ratios)
