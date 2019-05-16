CC=mpicc -g

simpar: simpar-mpi.o simpar.o
	$(CC) -o simpar simpar.o -lm
	$(CC) -o simpar-mpi simpar-mpi.o -lm

clean:
	rm *.o
	rm simpar-mpi
	rm simpar