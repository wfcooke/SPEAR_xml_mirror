#!/bin/csh
#------------------------------------------------------------------------------
#  netflx_ens.csh
#
#  DESCRIPTION: This is a script that is inteded to drive all the
#               pre-postprocessing stages of data manipulations on analysis nodes.
#               It is intended to be called by FRE at the "refineDiag" stage
#               which happens just before the components are post processed by frepp.
#               To make this happen a path to the (would be) script should appear
#               in the <refineDiag> tag of the xmls, e.g.,
#      <refineDiag script="$(NB_ROOT)/mom6/tools/analysis/MOM6_refineDiag.csh"/>
#               Note that the above script should exist when frepp is called.
#               This could be achieved by cloning the mom6 git repo in the <csh> section of the setup block
#               in the corresponding the gfdl platfrom. E.g.,
#        <csh><![CDATA[
#           source $MODULESHOME/init/csh
#           module use -a /home/John.Krasting/local/modulefiles
#           module purge
#           module load jpk-analysis/0.0.4
#           module load $(FRE_VERSION)
#         ]]></csh>
#
#------------------------------------------------------------------------------
echo ""
echo "  -- begin netflx.csh --  "
echo ""
#The mere existance of a refineDiag section in the xml pointing to any non-empty refineDiag
#script (as simple as doing a "echo" above)
#causes the history files to be unpacked by frepp in /ptmp/$USER/$ARCHIVE/$year.nc
#when the current year data lands on gfdl archive.
#
echo "  ---------- begin monthly analysis ----------  "
echo ""
#
#Generate this year's analysis figures based on the unpacked history files
#
#Try setting fre version to the caller version
if ( ! $?FREVERSION ) set FREVERSION = fre/bronx-19
set fremodule = $FREVERSION
set freanalysismodule = $FREVERSION

set src_dir=/home/wfc/SPEAR/SPEAR_include
# The following variables are set by frepp, but frepp is not called yet at refineDiag stage of FRE workflow,
#  so we need to explicitly set them here
set descriptor = $name
set out_dir = /home/wfc/SPEAR/SPEAR_include                  #Niki: How can we set this to frepp analysisdir /nbhome
set yr1 = $oname

# make sure valid platform and required modules are loaded
if (`gfdl_platform` == "hpcs-csc") then
   source $MODULESHOME/init/csh
   module use -a /home/fms/local/modulefiles /usr/local/paida/Modules
   module purge
   module load $fremodule
   module load $freanalysismodule
   module load gcc
   module load netcdf/4.2
#   module load python/2.7.3
#   module load python/2.7.12
else
   echo "ERROR: invalid platform"
   exit 1
endif

# check again?
if (! $?FRE_ANALYSIS_HOME) then
   echo "ERROR: environment variable FRE_ANALYSIS_HOME not set."
   exit 1
endif

#
#At this point of the FRE workflow we are in a /vftmp directory on an analysis node
#with all the history files for the current finished year ${yr1} unpacked and present
#
echo "We are inside the refineDiag script"
pwd
ls -l
# This script calculates the netflx at the surface from the formula
# netflx_sfc=swdn_sfc-swup_sfc+lwflx-2.5e6*evap-shflx

set script_dir=${out_dir}/refineDiag

set varlist=(`ls -1 *atmos_month.*.nc`)

foreach sne ( $varlist )
  # Get the tile number
  set tilename=$sne:r:e 
  # Check the variables are in the file
  ncks -m -v swup_sfc,swdn_sfc,lwflx,shflx,evap $sne > /dev/null
  if ($status == 0) then
    # Extract the 5 variables of interest to an interim file.
    ncks -v time_bnds,average_T1,average_T2,average_DT,swup_sfc,swdn_sfc,lwflx,evap,shflx $yr1.atmos_month.$tilename.nc interim.nc
    # Sum the 5 variables as needed. Output to a atmos_month_refined file.
    ncap2 -s "hfds=lwflx+swdn_sfc-swup_sfc-2.5e6*evap-shflx"  interim.nc interim1.nc 
    ncatted -h -a long_name,hfds,o,c,"net (down-up) heat flux at surface" interim1.nc
    ncatted -h -a _FillValue,hfds,o,d,1.0e+20 interim1.nc
    # Now remove the 5 variables from the refined file, leaving netflx_sfc,time_bnds and average_{T1,T2,DT}.
    ncks -x -v swup_sfc,swdn_sfc,lwflx,evap,shflx interim1.nc  $refineDiagDir/$yr1.atmos_month_refined.$tilename.nc
    rm -f interim.nc interim1.nc
  else
    echo "Missing one of lwflx, swdn_sfc, swup_sfc, evap, or shflx from input file."
    exit 1
  endif
end

echo "  ---------- end monthly analysis ----------  "

echo "  -- end netflx.csh --  "

exit 0
