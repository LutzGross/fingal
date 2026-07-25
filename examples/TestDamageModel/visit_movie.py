"""
VisIt CLI: render every load step of a cube run to a frame (mid-plane x2=L/2
damage slice) for a propagation movie. Frames go to <dir>/frames/frame%04d.png.

Run: /opt/visit-3.4.2/bin/visit -cli -nowin -s visit_movie.py [results_dir]
"""
import glob
import os
import sys

d = "cube_n24"
for a in sys.argv:
    if os.path.isdir(a):
        d = a
frames = os.path.join(d, "frames")
os.makedirs(frames, exist_ok=True)

OpenDatabase(os.path.join(d, "step*.silo database"))
AddPlot("Pseudocolor", "D")
p = PseudocolorAttributes()
p.minFlag, p.min = 1, 0.0
p.maxFlag, p.max = 1, 0.5
p.colorTableName = "hot"
SetPlotOptions(p)

AddOperator("Slice")
sa = SliceAttributes()
sa.originType = sa.Percent
sa.originPercent = 50
sa.axisType = sa.YAxis
sa.project2d = 1
SetOperatorOptions(sa)
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
for i in range(n):
    SetTimeSliderState(i)
    SaveWindow()

DeleteAllPlots()
quit()
