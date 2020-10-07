all: main

main: main.cpp poisson.cpp results.cpp
	g++ -O3 -std=c++11 -Wall -g3 -o $@ $^ -lpthread
	
%.prof: %
	-./$<
	gprof -b $< > $@

clean:
	rm main
	rm -f gmon.out perf.data*
