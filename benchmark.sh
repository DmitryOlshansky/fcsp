#!/bin/bash
TIMEFORMAT="%R";
for t in 1 2 3 4 5 6 7 8
do
        TIME=$(time (fcss-2a -t $t `cat list.txt` > /dev/null 2>&1) 2>&1)
        echo "$t,$TIME"
done
