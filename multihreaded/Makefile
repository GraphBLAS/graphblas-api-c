all:	example error

example:	example.cc graphblas.h
	g++ -g -fopenmp -o example -I. example.cc

error:	error.cc graphblas.h
	g++ -g -fopenmp -o error -I. error.cc

clean:
	rm example
