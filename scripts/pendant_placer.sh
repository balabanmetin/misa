#!/bin/bash

# $1 backbone tree
# $2 dist- col vector


s1=`awk '(NR==FNR){a[$1];next} ($1 in a)' <(nw_labels $1) $2 | sort -k2 -n | head -n 2 | cut -f1 |head -n 1`
s2=`awk '(NR==FNR){a[$1];next} ($1 in a)' <(nw_labels $1) $2  | sort -k2 -n | head -n 2 | cut -f1 |tail -n +2`


cat $1 | sed "s/$s1/($s1:0,mix_1:0)/g" | sed "s/$s2/($s2:0,mix_2:0)/g"



