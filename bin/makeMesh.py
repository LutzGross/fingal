#!/usr/bin/python3
"""
fingal - Creates a mesh fly file using the stations and - if available - a topography

"""
import argparse
import importlib, sys, os
from fingal import *
import numpy as np

sys.path.append(os.getcwd())



parser = argparse.ArgumentParser(description='Creates a mesh fly file using the station location information. the gmsh mesh generator is used.', epilog="fingal by l.gross@uq.edu.au version Feb/2024")#
parser.add_argument('--zlevel', '-z',  dest='zlevel', type=float, default=None, help="file name of the topography interpolation file.")
parser.add_argument('--core_depth', '-D',  dest='core_depth', type=float, default=40., help="core depth [%%] relative survey diameter.")
parser.add_argument('--extra_core', '-C',  dest='extra_core', type=float, default=40., help="extra core length [%%] relative to survey edge length.")
parser.add_argument('--extra_padding', '-P',  dest='extra_padding', type=float, default=300., help="extra padding length [%%] relative to core edge length.")
parser.add_argument('--num_elements_outer_edge', '-O',  dest='num_elements_outer_edge', type=int, default=10, help="numer of elements on outer edges.")
parser.add_argument('--mesh_size_core', '-S',  dest='mesh_size_core', type=float, default=None, help="element size in core region. If not given mean of distance to nearest electrode is used.")
parser.add_argument('--rel_mesh_size_electrodes', '-E',  dest='rel_mesh_size_electrodes', type=float, default=0.1, help="reduction factor of element size near electrodes relative to core mesh size.")
parser.add_argument('--topo', '-t',  dest='topofile', type=str, default=None, help="file name of the topography interpolation file.")
#parser.add_argument('--format', '-f', dest='topoformat', type=str, default = None, help = 'format of topography file ["npgrid", "npcloud" ]')
parser.add_argument('--plot', '-p',  dest='plot', type=str, default=None, help="file name to plot topography and stations in survey coordinate system.")
parser.add_argument('--recenter', '-R',  dest='recenter', action='store_true', default=False, help="recenter geometry to center of electrodes.")
parser.add_argument('--basename', '-B',  dest='basename', metavar="BASENAME", type=str, default="tmp", help="file names for intermediate files")
parser.add_argument('--1', '-1',  dest='step1', action='store_true', default=False, help="make flat surface geometry file BASENAME + _2D_flat.geo; then stop.")
parser.add_argument('--2', '-2',  dest='step2', action='store_true', default=False, help="make flat surface mesh file BASENAME + _2D_flat.msh; then stop.")
parser.add_argument('--3', '-3',  dest='step3', action='store_true', default=False, help="apply topography to surface mesh file BASENAME + _2D.msh; then stop.")
parser.add_argument('--4', '-4',  dest='step4', action='store_true', default=False, help="generate 3D geometry with topography BASENAME + .geo; then stop.")
parser.add_argument('--5', '-5',  dest='step5', action='store_true', default=False, help="generate 3D mesh file BASENAME + .msh; then stop.")
parser.add_argument('--6', '-6',  dest='step6', action='store_true', default=False, help="generate fly file FLYFILE.")
parser.add_argument('--fly', '-f',  dest='flyfile', metavar="FLYFILE", type=str, default=None, help="name of the fly file to be generated in coordinate system of the survey.")
parser.add_argument('--new_station_file', '-F',  dest='new_station_file', metavar="NEW_STATION_FILE", type=str, default=None, help="name of the new electrode file in coordinate system of the survey ('id, x,y, z').")
parser.add_argument(dest='stationfile', metavar='STATION_FILE', type=str, help='file of of electrodes with "id, x, y, z" or "id, x, y"')


args = parser.parse_args()

DELIMITER = ','
electrodes = {}
f = open(args.stationfile, 'r')
zlevel_e = 0.
zlevel_count =0
eDim=3
line=f.readline()
while line.strip():
        ll=line.strip().split(DELIMITER)
        S=int(ll[0])
        x=float(ll[1])
        y=float(ll[2])
        if len(ll) > 3 :
            z=float(ll[3])
            electrodes[S] = np.array([x, y, z])
            zlevel_e+=z
            zlevel_count+=1
        else:
            eDim=2
            electrodes[S]=np.array([x,y])
        line=f.readline()
f.close()
print(f"{len(electrodes)} read from {args.stationfile}.")

generator = MeshWithTopgraphy(electrodes=electrodes, recenter=args.recenter,
                              core_depth=args.core_depth,
                              extra_core=args.extra_core,
                              extra_padding=args.extra_padding,
                              num_elements_outer_edge=args.num_elements_outer_edge,
                              mesh_size_core=args.mesh_size_core,
                              rel_mesh_size_electrodes=args.rel_mesh_size_electrodes,
                            basename=args.basename)

if args.topofile:
    toptable= np.load(args.topofile)
    if 'x' in toptable.keys() and  'y' in toptable.keys()  and  'elevation' in toptable.keys():
        if toptable['elevation'].ndim == 2:
            generator.setTopgraphyFromGrid(x= toptable['x'], y= toptable['y'], elevation= toptable['elevation'] , method="nearest")
        elif toptable['elevation'].ndim == 1:
            generator.setTopgraphyFromCloud(x= toptable['x'], y= toptable['y'], elevation= toptable['elevation'] , method="nearest")
else:
    if args.zlevel:
        generator.setFlatTopography( zlevel = args.zlevel)
    else:
        if zlevel_count > 0 :
            generator.setFlatTopography( zlevel = zlevel_e/zlevel_count)
        else:
            raise ValueError("Unable to set z-level. Use option --zlevel ")

if args.plot:
    generator.plotting(args.plot)
else:
    run_all = not any( [args.step1, args.step2, args.step3, args.step4, args.step5, args.step6] )
    if  (args.step6 or run_all) and not args.flyfile:
        raise ValueError("No fly file given but generation request. use --fly option.")

    if args.step1 or run_all:
        generator.makeFlatSurfaceGeometry()
    if args.step2 or run_all:
        generator.makeFlatSurfaceMesh()
    if args.step3 or run_all:
        generator.addTopologyTo2DMesh()
    if args.step4 or run_all:
        generator.generate3DGeometry()
    if args.step5 or run_all:
        generator.generate3DMesh()
    if (args.step6 or run_all) and args.flyfile:
        generator.toFlyFile(flyfile=args.flyfile)
if args.new_station_file:
    generator.writeStations(args.new_station_file, delimiter = DELIMITER)
print("All done - Have a nice day!")



