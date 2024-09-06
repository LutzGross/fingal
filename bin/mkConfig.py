#!/usr/bin/python3
"""
fingal - program to set an initial configuration file.
"""
import argparse, os
from datetime import datetime
from fingal import SurveyData

parser = argparse.ArgumentParser(description='This creates an initial configuration file that needs further editing.')
parser.add_argument(dest='project', metavar='PROJECT', type=str, help='configuration file name (py extension will be added.)')
parser.add_argument('--station', '-s',  dest='stationfile', type=str, default=None, help="name of station file")
parser.add_argument('--data', '-d',  dest='datafile', type=str, default=None, help="name of data file")
parser.add_argument('--ipole', '-i',  dest='ipole', action='store_true', default=False, help="injections are pole")
parser.add_argument('--mpole', '-m',  dest='mpole', action='store_true', default=False, help="measurements are pole")
parser.add_argument('--obs', '-o',  dest='obs', type=str, default='R,ETA', help="descriptor of columns in data file ")


args = parser.parse_args()

out=""
for line in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "config_template.txt"), 'r').readlines():
        if not (args.compact and line.startswith("#.")):
            out+=line
datacolumns=[]
for o in args.obs.split(','):
    assert SurveyData.checkObservationType(o), 'Unknown observation type {o}'
    datacolumns.append(o)
    
CHANGE= {
"created" :  datetime.now().strftime("%d.%m.%Y %H:%M"),
"project" : args.project,
"meshfile" : args.project+".fly",
"stationfile" : args.stationfile,
"stationsFMT" : "s%s",
"datafile" : args.datafile,
"dipoleInjections" : not args.ipole,
"dipoleMeasurements" : not args.mpole,
"datacolumns" :datacolumns,
"schedulefile" : None,
"core" : ['core'],
"face" : ['faces']
}


FN=args.project+".py"
open(FN, 'w').write(out.format(**CHANGE))
print(f"Configuration file {FN} created.")
print(f"Please check the content of the file.")

