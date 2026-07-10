export OMP_PROC_BIND=spread
export OMP_PLACES=threads

for p in 1 2 3 4 5 6 7 8
do
  export OMP_NUM_THREADS=$p
  ./minetest.py  2>&1 > "logs/log$p"
done