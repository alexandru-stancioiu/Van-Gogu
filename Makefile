IN="in/julia2.in"
OUT="julia2.out"

build:
	mpicc -o fractal fractal.c -lm
run:
	mpirun -np 4 ./fractal $(IN) $(OUT)
clean:
	rm -rf *~ fractal *.out
