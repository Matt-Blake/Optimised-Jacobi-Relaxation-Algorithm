all: main

main: main.cpp poisson.cpp results.cpp
	g++ -O3 -std=c++11 -Wall -g3 -pg -o $@ $^ -lpthread
	
%.prof: %
	-./$< 
	gprof -b $< > $@

clean:
	rm -f main
	rm -f gmon.out gmon.sum main perf.data*
