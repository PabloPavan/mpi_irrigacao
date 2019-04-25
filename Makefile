# comentario -> para executar digite make 

all: compile
	
compile:
	mpicc -o irrigacao-mpi irrigacao-mpi.c


clean:
	rm -rf *.out

