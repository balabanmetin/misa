#!/bin/bash

# $1 placement tree
# $2 things.txt
# $3 true.tree

s1=`head -n 1 $2`
s2=`tail -n +2 $2`

t1=`cat $1 | sed "s/mix_1/$s1/g" `
t2=`cat $1 | sed "s/mix_2/$s2/g" `
t3=`cat $1 | sed "s/mix_1/$s2/g" `
t4=`cat $1 | sed "s/mix_2/$s1/g" `

v1=`compareTrees.missingBranch $3 <(echo $t1) -simplify | cut -f2 -d' '`
v2=`compareTrees.missingBranch $3 <(echo $t2) -simplify | cut -f2 -d' '`
v3=`compareTrees.missingBranch $3 <(echo $t3) -simplify | cut -f2 -d' '`
v4=`compareTrees.missingBranch $3 <(echo $t4) -simplify | cut -f2 -d' '`



if ((v1+v2 < v3+v4 )); then
    printf "%s\n%s\n" $v1 $v2
else
    printf "%s\n%s\n" $v3 $v4
fi
