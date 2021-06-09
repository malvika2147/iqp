#!/bin/bash

for n in {40..300..10};
do
    echo $n >in;
    echo ${n};
    filename=data${n};
    echo -n >$filename;
    time for i in {1..200}; do ./main <in >>$filename; done
done

