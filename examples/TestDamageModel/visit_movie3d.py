"""
VisIt CLI: 3-D propagation movie of a cube run. A Threshold operator hides the
undamaged bulk (D < DMIN) and the surviving damaged region is coloured by D
(hot, [0, 0.5]); the camera orbits so the growing 3-D band is legible. Frames go
to <dir>/frames3d/frame%04d.png.

Run: /opt/visit-3.4.2/bin/visit -cli -nowin -s visit_movie3d.py [results_dir]
"""
import glob
import math
import os
import sys

DMIN = 0.1                                       # hide material below this damage

d = "cube_n24"
for a in sys.argv:
    if os.path.isdir(a):
        d = a
frames = os.path.join(d, "frames3d")
os.makedirs(frames, exist_ok=True)

OpenDatabase(os.path.join(d, "step*.silo database"))
AddPlot("Pseudocolor", "D")
p = PseudocolorAttributes()
p.minFlag, p.min = 1, 0.0
p.maxFlag, p.max = 1, 0.5
p.colorTableName = "hot"
SetPlotOptions(p)

AddOperator("Threshold")
t = ThresholdAttributes()
t.listedVarNames = ("D",)
t.lowerBounds = (DMIN,)
t.upperBounds = (1.0e37,)
SetOperatorOptions(t)
DrawPlots()
ResetView()

sw = SaveWindowAttributes()
sw.family = 1
sw.outputToCurrentDirectory = 0
sw.outputDirectory = frames
sw.fileName = "frame"
sw.format = sw.PNG
sw.width, sw.height = 800, 800
SetSaveWindowAttributes(sw)

n = TimeSliderGetNStates()
print("timesteps:", n)
el = math.radians(22.0)               # camera elevation
a0 = math.radians(-125.0)             # start azimuth
sweep = math.radians(80.0)            # total orbit over the movie
for i in range(n):
    SetTimeSliderState(i)
    a = a0 + sweep * i / max(n - 1, 1)
    v = GetView3D()
    v.viewNormal = (math.cos(el) * math.cos(a),
                    math.cos(el) * math.sin(a),
                    math.sin(el))
    v.viewUp = (0.0, 0.0, 1.0)
    SetView3D(v)
    SaveWindow()

DeleteAllPlots()
quit()
