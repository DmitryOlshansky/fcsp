#!/bin/bash

for t in tests/* ; do
	echo "Encoding" `echo -n $t | sed -r 's|.*/(.*)|\1|'`
	FILES=`echo $t/MOL/*.{mol,MOL} | sort`
	./fcss-2a --format=csv -d . $FILES > $t/fcss-2a-dev.csv
done 2>test-suite.log

# Sample of test out
head -6 test-suite.log
echo ....
if [ `wc -l test-suite.log | cut -f1 -d ' '` -gt 6 ] ; then
	tail -6 test-suite.log
fi
