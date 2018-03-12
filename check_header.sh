#!/bin/bash
ls -l *.h > hlist_new
RET=0
if diff hlist hlist_new >> /dev/null; then
    RET=1
fi
mv hlist_new hlist -f
echo $RET
