# comentario -> para executar digite make 

all: compile run
	
compile:
	mpicc irrigacao-mpi.c -fopenmp -lm

run:
	mpirun --mca btl ^openib --mca btl_tcp_if_include enp7s0f1 --bind-to none -np 88 ./a.out

clean:
	rm -rf a.out

