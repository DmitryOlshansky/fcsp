#!/bin/bash

for t in tests/* ; do
	echo "Encoding" `echo -n $t | sed -r 's|.*/(.*)|\1|'`
	FILES=`echo $t/MOL/*.{mol,MOL} | sort`
	./fcss-2a --format=csv -d . $FILES > $t/fcss-2a-dev.csv
done 2>test-suite.log

# Sample of test out
if [ `wc -l test-suite.log | cut -f1 -d ' '` -gt 12 ] ; then
	head -6 test-suite.log
	echo ....
	tail -6 test-suite.log
else
	cat test-suite.log
fi

for t in tests/* ; do
	echo "Comparing " `echo -n $t | sed -r 's|.*/(.*)|\1|'`
	./fcss-comp $t/fcss-2.csv $t/fcss-2a-dev.csv | grep -A 5 "SUMMARY" | tail -2
done
