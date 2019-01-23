#! /bin/bash --norc
#help and version handling
ixit=0
for i in $*; do tst=$i #echo $tst
if [ "${tst#-}" != $tst ] ; then ixit=1
  if [ "${tst##--}" == version ] ;  then
    echo "Introduces/interacts with SMARDDA-LIB 1.0"
    echo "LGPL 3 2019"
  elif [ "${tst##--}" == help ] ;  then
    echo "Introduces the occasional linux user to ray tracing"
    echo "using SMARDDA-LIB"
    echo "with web-links to extensive more advanced material."
  elif [ "${tst##--}" == tmp ] ;  then
    ixit=2
  elif [ "${tst##--}" == force ] ;  then
    ixit=3
  else
    echo "Valid options are --help and --version"
    echo "EXPERT options are --tmp and --force"
  fi
fi
done
if [ $ixit == 2 ] ; then exit 2;fi
if [ $ixit == 3 ] ; then exit 3;fi
if [ $ixit != 0 ] ; then exit 1;fi
