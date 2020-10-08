#!/bin/bash


echo "Bash version ${BASH_VERSION}..."
echo "Terminal output save test"

make clean
make all
break>thread_output.txt

for x in 101 501 901
do
	for i in {1..10}
	do
		for j in {1..10}
			do
				./main $x 10 $i >> thread_output.txt
			done
	done
done

