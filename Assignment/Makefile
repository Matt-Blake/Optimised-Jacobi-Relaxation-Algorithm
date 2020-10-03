all: poisson_test

poisson_test: poisson_test.cpp poisson.cpp
	g++ -O3 -std=c++11 -Wall -g3 -o $@ $^

clean:
	rm poisson_test
