#!/usr/bin/python3

out = "export OMP_NUM_THREADS=6\n"

jj=0
jj0=jj
for contrast in [1., 100.]:
    for k in [1, 2, 4, 8, 16, 32, 64, 128, 256]:
        for ratio in [0., 0.02, 0.08, 0.32, 0.64, 1.28 ]:
            line = f"./unstructest.py  --contrast {contrast} --ratio {ratio}  meshes/dom{k}.fly 2>&1 > logs/log{jj}"
            out+= line+"\n"
            jj+=1
with open("runs.sh","w") as f:
    f.write(out)

print(f"{jj-jj0} tests created.")