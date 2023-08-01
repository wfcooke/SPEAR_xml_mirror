Repository for SPEAR model.

This contains the XML and related files to run the SPEAR model.

module load fre/bronx-20
git clone -b 2023.02 https://gitlab.gfdl.noaa.gov/SPEAR/xml.git
multi-fremake -x SPEAR_c192_o1_Historical_IC1921_Q50_ensembles.xml -p ncrc5.intel22 -t repro-openmp,prod-openmp -e SPEAR_Q2023.02_nonsymMOM6_exec

! Experiments:
  1. SPEAR_c192_o1_Hist_AllForc_IC1921_Q50_ens_01_03 will give you a test of a 3 member ensemble.
  2. SPEAR_c192_o1_Control_1850_Q50.xml will give you a test of a Control (fixed 1850 forcing) run.

multi-frerun -x SPEAR_c192_o1_Historical_IC1921_Q50_ensembles.xml -p ncrc5.intel22 -t prod-openmp,repro-openmp -r basic_prodsettings --overwrite --no-transfer -s -e SPEAR_c192_o1_Hist_AllForc_IC1921_Q50_ens_01_03
multi-frerun -x SPEAR_c192_o1_Control_1850_Q50.xml -p ncrc5.intel22 -t repro-openmp,prod-openmp -r basic_prodsettings --overwrite -s -e SPEAR_c192_o1_Control_1850_Q50
