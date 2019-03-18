#! /bin/bash
# testdeck for smardda-lib, assumes software installed under directory HS
# The test outputs will go in directory $HS/TEST
#
if [[ -n "$SMALIB_DIR" ]] ; then HS=$SMALIB_DIR; else HS=$HOME/smardda/smardda-lib; fi
export HS
wd=$PWD
#
# geofil test
#copy $SMALIB_DIR/TEST/GFRUN to a work
#directory, removing the comparison output files, run geofil, then 
#compare results, ie:
cd $HS/TEST
if [[ -d workdirg ]] ; then rm -rf workdirg ; fi
cp -r GFRUN workdirg
cd workdirg
rm -f absb_geofil.out absb.vtk absb.log
geofil absb
diff -b absb.vtk ../GFRUN/absb.vtk
# visualise results using ParaView
#
# generate tracks
#To produce a set of 400 points defining 200 tracks as adjacent vectors
#in the file exps00000200.qry
cd $HS/TEST
if [[ -d workdirr ]] ; then rm -rf workdirr ; fi
cp -r RAYS workdirr
cd workdirr
if [[ ! -e exps ]] ; then
   if [[ $(command -v ifort) ]] ; then FORT=ifort ; else FORT=gfortran ; fi
   if [[ -e ../../fortd/libg.a ]] ; then
   $FORT -g exps.f -o exps ../../fortd/libg.a ; else 
   $FORT exps.f -o exps ../../fortd/lib.a ; fi
fi
./exps 200 4 0
diff -b exps00000200.qry ../RAYS/exps00000200.qry
# hdsgen and move test
#copy $SMALIB_DIR/TEST/MTRUN to a work
#directory, removing the comparison output files, run hdsgen and move, then 
#compare results, ie:
cd $HS/TEST
if [[ -d workdirm ]] ; then rm -rf workdirm ; fi
cp -r MTRUN workdirm
cd workdirm
rm -f *_* *.log
hdsgen testhds
diff -b testhds_hds.hds ../MTRUN/testhds_hds.hds
move test
diff -b test_movx.vtk ../MTRUN/test_movx.vtk
# visualise results using ParaView
cd $wd
