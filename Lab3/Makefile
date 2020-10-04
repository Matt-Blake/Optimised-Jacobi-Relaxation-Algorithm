CFLAGS=-g -Wall -pipe -lm

# Turn on the gperf profiling
CFLAGS+=-pg

CC=g++

profile_example: profile_example.cpp
	$(CC) -o $@ $< $(CFLAGS) $(LFLAGS)

%.prof: %
	-./$<
	gprof -b $< > $@

clean:
	rm -f profile_example
	rm -f gmon.out perf.data*
