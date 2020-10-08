#!/bin/bash


echo "Bash version ${BASH_VERSION}..."
echo "Terminal output save test"

make clean
make all
break>thread_output.txt

for x in 101 201 301 401 501 601 701 801 901
do
	for i in {1..8}
	do
		for j in {1..10}
			do
				./main $x 10 $i >> thread_output.txt
			done
	done
done
