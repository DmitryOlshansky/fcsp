#!/bin/bash

echo "Выборка,Только в ФКСП-2,Только в ФКСП-2А,Общие для ФКСП-2 и ФКСП-2А"
for t in tests/* ; do
	DATA=`fcss-comp -i $t/fcss-2.csv $t/fcss-2a.csv |  grep -A 5 "SUMMARY" | tail -3`
	COUNTS=`echo "$DATA" | cut -f2 -d ':'`
	echo -n $t
	for f in $COUNTS ; do
		echo -n ",$f"
	done
	echo
done
