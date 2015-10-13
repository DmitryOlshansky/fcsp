#!/bin/bash

echo "file,missing,extra"
for t in tests/* ; do
	DATA=`fcss-comp $t/fcss-2.csv $t/fcss-2a.csv |  grep -A 5 "SUMMARY" | tail -2`
	MISSING=`echo "$DATA" | sed -r '/EXTRA:.*/d;s/MISSING:\s*(.*)/\1/'`
	EXTRA=`echo "$DATA" | sed -r '/MISSING:.*/d;s/EXTRA:\s*(.*)/\1/'`
	echo $t, $MISSING, $EXTRA
done
