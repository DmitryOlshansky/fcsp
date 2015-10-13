#!/bin/bash
LOG_LEVEL=3 # warn-s and worse
for t in tests/* ; do
	echo "Encoding" `echo -n $t | sed -r 's|.*/(.*)|\1|'`
	FILES=`echo $t/MOL/*.MOL | sort`
#	echo $FILES
	echo "$FILES" | xargs -n 2 ./fcss-2a -v ${LOG_LEVEL} --format=csv  > $t/fcss-2a-dev.csv
done 2>test-suite.log

for t in tests/* ; do
	echo "Comparing " `echo -n $t | sed -r 's|.*/(.*)|\1|'`
	./fcss-comp -i $t/fcss-2.csv $t/fcss-2a-dev.csv | tee $t/diff.cmp | grep -A 5 "SUMMARY" | tail -3
done
