#!/bin/bash

# $1 backbone tree
# $2 dist mat

source activate three;

tmp=`mktemp -t "XXXXXX.jplace"`
tmp2=`mktemp -t "XXXXXX.nwk"`

python3 ~/apples/run_apples.py -t $1  -d $2 -o $tmp -T 1
guppy tog -o $tmp2 $tmp
cat $tmp2 | sed "s/mix:/(mix_1:0,mix_2:0):/g"
rm $tmp $tmp2




