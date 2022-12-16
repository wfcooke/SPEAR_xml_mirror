#!/bin/csh

set xmlDir = $rtsxml:h
set refineDiagScriptDir = $xmlDir/SPEAR_include/refineDiag 
set splitter=$refineDiagScriptDir/comp_ensemble_avg_v4_fre_atmos_only

set startDate = `date`
echo $startDate
set startDateS = `date '+%s'`

echo Using environment:
env
set

# Set some defaults
if ( ${?platform} == 0 ) then
   setenv FRE_PLATFORM default
else
   setenv FRE_PLATFORM ${platform}
endif
if ( ${?target} == 0 ) then
   setenv FRE_TARGET default
else
   setenv FRE_TARGET ${target}
endif

echo Beginning ensemble average calculation and ensemble member frepp submission
echo $startDate

$splitter -v -y .true. \
          -b $oname \
          -c .true. \
          -x $rtsxml \
          -s .true. \
          -p .true. \
          -d $archive/history \
          -o $archive/ensemble \
          -w $work \
          -P $ptmpDir \
          $name

set errorCode=$?
set endDate=`date`
set endDateS=`date '+%s'`
@ time = $endDateS - $startDateS

if ( $errorCode != 0 ) then
   echo "Error while running $splitter, status $errorCode"
   echo $endDate
   exit $errorCode
endif

echo Completed ensemble average calculation and ensemble member frepp submission in $time seconds
echo $endDate

echo "Workdir is $work"
ls $work

# Frepp will stop all processing of the refineDiag script if the exit
# value is less than 0.
exit -1
