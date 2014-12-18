#!/bin/csh

if ($#argv < 4) then
  echo "Usage: getdata.csh list_file gapO gapE bias"
  echo " V.1.0.0 by Yao-ming Huang 2005/8/9"
  exit
endif

if ($#argv >= 4) then
  setenv LISTFILE $1
  setenv GAPO $2
  setenv GAPE $3
  setenv BIAS $4
endif

setenv RUNDIR `pwd`

@ N = 0
@ MX = `wc -l $RUNDIR/$LISTFILE | awk '{print $1}'`
while ($N < $MX)
  @ N ++
  setenv FILE1 `head -$N $RUNDIR/$LISTFILE | tail -1 | awk '{print $1}'`
  setenv FILE2 `head -$N $RUNDIR/$LISTFILE | tail -1 | awk '{print $2}'`
  cd $RUNDIR/gca
  $RUNDIR/gca/xalignD $FILE1 $FILE2 $GAPO $GAPO $GAPE $GAPE $BIAS $BIAS
  sleep 1
end
while ( `cat $RUNDIR/gca/*exp.txt | grep "Aver" | wc -l | awk '{print $1}'` < $MX)
  sleep 5
end
rm $RUNDIR/ref/*_exp.txt
mv $RUNDIR/gca/*_exp.txt $RUNDIR/ref
