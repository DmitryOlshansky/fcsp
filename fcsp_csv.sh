#!/bin/sh
NAME=/tmp/fcsp.codes
DESC_DIR="." # load descriptors from this directory
if [ $# -lt  1 ] ; then
	exit 1
fi
for f in $@ ; do
	MOL=`echo $f | sed -r 's|.*/(.*)|\1|'`
	./fcsp -d $DESC_DIR $f > $NAME 2>/tmp/fcsp.log
	FILE=`echo $f | sed -r 's|.*/(.*)|\1|'`
	(echo $MOL; wc -l $NAME | sed -r 's/([0-9]+)\s*.*/\1/'; (cat $NAME | tr '\n' ' ')) | tr '\n' ';'
	echo
done
