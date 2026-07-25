"""
VisIt CLI: render a chosen scalar field from a cube run's last silo -- a 3-D
surface view and a vertical mid-plane (x2=L/2) slice, colour range autoscaled.

Run: /opt/visit-3.4.2/bin/visit -cli -nowin -s visit_field.py <results_dir> <var>
     (var defaults to "D"; e.g. "phi" for porosity, "sig_zz" for axial stress)
"""
import glob
import os
import sys

d, var = "cube_n24", "D"
for a in sys.argv:
    if os.path.isdir(a):
        d = a
    elif not a.startswith("-") and not a.endswith(".py"):
        var = a

files = sorted(glob.glob(os.path.join(d, "step*.silo")))
last = files[-1]
OpenDatabase(last)
md = GetMetaData(last)
print("scalars:", [md.GetScalars(i).name for i in range(md.GetNumScalars())])
print("rendering", var, "from", last)

AddPlot("Pseudocolor", var)
p = PseudocolorAttributes()
p.minFlag, p.maxFlag = 0, 0            # autoscale to the data range
p.colorTableName = "hot"
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
sw.fileName = "visit_%s_3d" % var
SetSaveWindowAttributes(sw)
SaveWindow()

AddOperator("Slice")
sa = SliceAttributes()
sa.originType = sa.Percent
sa.originPercent = 50
sa.axisType = sa.YAxis
sa.project2d = 1
SetOperatorOptions(sa)
DrawPlots()
ResetView()

sw.fileName = "visit_%s_slice" % var
SetSaveWindowAttributes(sw)
SaveWindow()

DeleteAllPlots()
CloseDatabase(last)
quit()
