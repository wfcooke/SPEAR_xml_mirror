Repository for SPEAR model.

This contains the XML and related files to run the SPEAR model.

```git
git clone -b 2022.01_pruned https://gitlab.gfdl.noaa.gov/SPEAR/xml.git
```

```
fremake -x SPEAR_c192_o1_Historical_IC1921_Q50_ensembles.xml -p ncrc4.intel18 -t prod-openmp SPEAR_Q2022.01_nonsymMOM6_exec
```

! Run experiments:
  1. SPEAR_c96_o1_Control_1850_Q50 will give you a test of a Control (fixed 1850 forcing) run.
  2. SPEAR_c192_o1_Control_1850_Q50 will give you a test of a Control (fixed 1850 forcing) run.
  3. SPEAR_c384_o1_Control_1850_Q66_IC101_1PctTo2X will give you a test of a Control (fixed 1850 forcing) run.
  4. SPEAR_c192_o1_Hist_AllForc_IC1921_Q50_ens_01_03 will give you a test of a 3 member ensemble.

```
frerun -x SPEAR_c96_o1_Control_1850_Q50.xml -p ncrc4.intel18 -t prod-openmp -r basic SPEAR_c96_o1_Control_1850_Q50 --no-transfer -s
frerun -x SPEAR_c192_o1_Historical_IC1921_Q50_ensembles.xml -p ncrc4.intel18 -t prod-openmp -r basic_prodsettings SPEAR_c192_o1_Hist_AllForc_IC1921_Q50_ens_01_03 --overwrite --no-transfer -s
frerun -x SPEAR_c192_o1_Control_1850_Q50.xml -p ncrc4.intel18 -t prod-openmp -r basic_prodsettings SPEAR_c192_o1_Control_1850_Q50 --overwrite --no-transfer -s
frerun -x SPEAR_c384_o1_Control_1850_Q66_IC101_1PctTo2X.xml -p ncrc4.intel18 -t prod-openmp -r basic SPEAR_c384_o1_Control_1850_Q66_IC101_1PctTo2X --no-transfer -s
```
