#!/bin/sh
NAME=/tmp/fcsp.codes
./fcsp -d . $1 > $NAME 2>/tmp/fcsp.log
FILE=`echo $1 | sed -r 's|.*/(.*)|\1|'`
(echo $FILE; wc -l $NAME | sed -r 's/([0-9]+)\s*.*/\1/'; (cat $NAME | tr '\n' ' ')) | tr '\n' ';'
