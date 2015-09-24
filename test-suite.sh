#!/bin/bash
for t in tests/* ; do
	echo "Encoding" `echo -n $t | sed -r 's|.*/(.*)|\1|'`
	./fcsp --format=csv -d . $t/MOL/*.{mol,MOL} > $t/fccs-2a-dev.csv
done 2>test-suite.log
