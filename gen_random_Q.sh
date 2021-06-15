#!/bin/bash

n=1000
echo $n $n >in;
filename=data${n}x${n};
echo -n >$filename;
for i in {1..200}; do 
    echo $i;
    time ./main <in >>$filename; 
done


