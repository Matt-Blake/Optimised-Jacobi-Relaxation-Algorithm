#!/bin/bash

echo "Profiling Commenced"


make clean
make all

./main 101 10 4
mv gmon.out gmon.sum
./main 101 10 4
gprof -s main gmon.out gmon.sum
./main 101 10 4
gprof -s main gmon.out gmon.sum
./main 101 10 4
gprof -s main gmon.out gmon.sum
./main 101 10 4
gprof -s main gmon.out gmon.sum
gprof main gmon.sum > profile-101.txt
echo "Created profile-101.txt"


make clean
make all

./main 201 10 4
mv gmon.out gmon.sum
./main 201 10 4
gprof -s main gmon.out gmon.sum
./main 201 10 4
gprof -s main gmon.out gmon.sum
./main 201 10 4
gprof -s main gmon.out gmon.sum
./main 201 10 4
gprof -s main gmon.out gmon.sum
gprof main gmon.sum > profile-201.txt
echo "Created profile-201.txt"


make clean
make all

./main 301 10 4
mv gmon.out gmon.sum
./main 301 10 4
gprof -s main gmon.out gmon.sum
./main 301 10 4
gprof -s main gmon.out gmon.sum
./main 301 10 4
gprof -s main gmon.out gmon.sum
./main 301 10 4
gprof -s main gmon.out gmon.sum
gprof main gmon.sum > profile-301.txt
echo "Created profile-301.txt"


make clean
make all

./main 401 10 4
mv gmon.out gmon.sum
./main 401 10 4
gprof -s main gmon.out gmon.sum
./main 401 10 4
gprof -s main gmon.out gmon.sum
./main 401 10 4
gprof -s main gmon.out gmon.sum
./main 401 10 4
gprof -s main gmon.out gmon.sum
gprof main gmon.sum > profile-401.txt
echo "Created profile-401.txt"


make clean
make all

./main 501 10 4
mv gmon.out gmon.sum
./main 501 10 4
gprof -s main gmon.out gmon.sum
./main 501 10 4
gprof -s main gmon.out gmon.sum
./main 501 10 4
gprof -s main gmon.out gmon.sum
./main 501 10 4
gprof -s main gmon.out gmon.sum
gprof main gmon.sum > profile-501.txt
echo "Created profile-501.txt"


make clean
make all

./main 601 10 4
mv gmon.out gmon.sum
./main 601 10 4
gprof -s main gmon.out gmon.sum
./main 601 10 4
gprof -s main gmon.out gmon.sum
./main 601 10 4
gprof -s main gmon.out gmon.sum
./main 601 10 4
gprof -s main gmon.out gmon.sum
gprof main gmon.sum > profile-601.txt
echo "Created profile-601.txt"


make clean
make all

./main 701 10 4
mv gmon.out gmon.sum
./main 701 10 4
gprof -s main gmon.out gmon.sum
./main 701 10 4
gprof -s main gmon.out gmon.sum
./main 701 10 4
gprof -s main gmon.out gmon.sum
./main 701 10 4
gprof -s main gmon.out gmon.sum
gprof main gmon.sum > profile-701.txt
echo "Created profile-701.txt"


make clean
make all

./main 801 10 4
mv gmon.out gmon.sum
./main 801 10 4
gprof -s main gmon.out gmon.sum
./main 801 10 4
gprof -s main gmon.out gmon.sum
./main 801 10 4
gprof -s main gmon.out gmon.sum
./main 801 10 4
gprof -s main gmon.out gmon.sum
gprof main gmon.sum > profile-801.txt
echo "Created profile-801.txt"


make clean
make all

./main 901 10 4
mv gmon.out gmon.sum
./main 901 10 4
gprof -s main gmon.out gmon.sum
./main 901 10 4
gprof -s main gmon.out gmon.sum
./main 901 10 4
gprof -s main gmon.out gmon.sum
./main 901 10 4
gprof -s main gmon.out gmon.sum
gprof main gmon.sum > profile-901.txt
echo "Created profile-901.txt"



