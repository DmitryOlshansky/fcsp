#!/bin/bash

for t in tests/* ; do
	echo "Encoding" `echo -n $t | sed -r 's|.*/(.*)|\1|'`
	FILES=`echo $t/MOL/*.{mol,MOL} | sort`
	./fcss-2a --format=csv -d . $FILES > $t/fcss-2a-dev.csv
done 2>test-suite.log

for t in tests/* ; do
	echo "Comparing " `echo -n $t | sed -r 's|.*/(.*)|\1|'`
	./fcss-comp $t/fcss-2.csv $t/fcss-2a-dev.csv | grep -A 5 "SUMMARY" | tail -2
done
