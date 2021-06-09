#!/bin/bash

n=50
for m in {10..100..5};
do
    echo $m $n >in;
    echo ${m}x${n};
    filename=data${m}x${n};
    echo -n >$filename;
    time for i in {1..100}; do ./main <in >>$filename; done
done

