#!/bin/bash

echo "Выборка,Соединений,Общие Дескрипторы,Только в ФКСП-2А,Только в ФКСП-2"
for t in tests/* ; do
	MOLS=`find $t/MOL -name '*.MOL' | wc -l`
	DATA=`fcss-comp -i $t/fcss-2.csv $t/fcss-2a.csv |  grep -A 5 "SUMMARY" | tail -3`
	COUNTS=`echo "$DATA" | cut -f2 -d ':'`
	echo -n $t,$MOLS
	for f in $COUNTS ; do
		echo -n ",$f"
	done
	echo
done
