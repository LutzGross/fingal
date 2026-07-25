"""
VisIt CLI script: render the damage field of a cube run from its silo files.
Produces a 3-D surface view (the shear-band X-pattern, cf. Fig. 7) and a vertical
mid-plane slice. Run with:

    /opt/visit-3.4.2/bin/visit -cli -nowin -s visit_damage.py [results_dir]
"""
import glob
import os
import sys

d = "cube_n24"
for a in sys.argv:
    if os.path.isdir(a):
        d = a
files = sorted(glob.glob(os.path.join(d, "step*.silo")))
last = files[-1]
print("rendering", last)

OpenDatabase(last)
md = GetMetaData(last)
print("scalars:", [md.GetScalars(i).name for i in range(md.GetNumScalars())])

AddPlot("Pseudocolor", "D")
p = PseudocolorAttributes()
p.minFlag, p.min = 1, 0.0
p.maxFlag, p.max = 1, 1.0
p.colorTableName = "rainbow"
SetPlotOptions(p)
DrawPlots()

v = GetView3D()
v.viewNormal = (-0.6, -0.7, 0.4)
v.viewUp = (0.0, 0.0, 1.0)
SetView3D(v)

sw = SaveWindowAttributes()
sw.family = 0
sw.outputToCurrentDirectory = 0
sw.outputDirectory = d
sw.format = sw.PNG
sw.width, sw.height = 900, 800
sw.fileName = "visit_damage_3d"
SetSaveWindowAttributes(sw)
SaveWindow()

# vertical mid-plane slice (x2 = L/2 -> x1-x3 plane through the weak zone)
AddOperator("Slice")
sa = SliceAttributes()
sa.originType = sa.Percent
sa.originPercent = 50
sa.axisType = sa.YAxis
sa.project2d = 1
SetOperatorOptions(sa)
DrawPlots()
ResetView()

sw.fileName = "visit_damage_slice"
SetSaveWindowAttributes(sw)
SaveWindow()

DeleteAllPlots()
CloseDatabase(last)
quit()                       # clean VisIt CLI exit (sys.exit segfaults)
