#!/bin/bash


echo "Bash version ${BASH_VERSION}..."
echo "Terminal output save test"

make clean
make all
break>compiler_output.txt

for j in {1..30}
do
	./main 501 10 4 >> compiler_output.txt
done

