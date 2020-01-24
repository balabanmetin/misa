#!/bin/bash

# $1 placement tree
# $2 things.txt
# $3 true.tree

s1=`head -n 1 $2`
s2=`tail -n +2 $2`

t1=`cat $1 | sed "s/mix_1/$s1/g" | sed "s/mix_2/$s2/g"`
t2=`cat $1 | sed "s/mix_2/$s1/g" | sed "s/mix_1/$s2/g"`

v1=`compareTrees.missingBranch $3 <(echo $t1) -simplify | cut -f2 -d' '`
v2=`compareTrees.missingBranch $3 <(echo $t2) -simplify | cut -f2 -d' '`

printf "%s\n%s\n" $v1 $v2 | sort -n | head -n 1


