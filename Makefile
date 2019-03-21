CC=gcc

simpar: simpar-omp.o simpar.o
	$(CC) -o simpar simpar.o -lm
	$(CC) -o simpar-omp simpar-omp.o -lm

clean:
	rm *.o
	rm simpar-omp
	rm simpar