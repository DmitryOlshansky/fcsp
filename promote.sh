#!/bin/bash
# Promote current dev encoding results as the new standard results
# Makes sense after each improvement and confirmed bugfix
for t in tests/* ; do
	cp $t/fcss-2a-dev.csv $t/fcss-2a.csv
done
echo "Done. See git status:" >&2
git status