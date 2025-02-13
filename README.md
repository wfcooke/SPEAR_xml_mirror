## Repository for SPEAR model

This contains the XML and related files to run the regression testing for the SPEAR model.

The legacy xml/experiments use the ascii diag, data, and field tables

The canopy experiments use the fre/2025.01 to build "bare metal" executables and containers

```
module load fre/bronx-23
git clone -b RTS https://gitlab.gfdl.noaa.gov/SPEAR/xml.git
multi-fremake -x SPEAR.xml -p ncrc5.intel23-classic -t repro-openmp,prod-openmp -e SPEAR_Q2025.01_nonsymMOM6_exec
```
Experiments:
  1. SPEAR_c192_o1_Hist_AllForc_IC1921_Q50_ens_01_03 will give you a test of a 3 member ensemble.
  2. SPEAR_c192_o1_Control_1850_Q50.xml will give you a test of a Control (fixed 1850 forcing) run.

```
multi-frerun -x SPEAR.xml -p ncrc5.intel23-classic -t prod-openmp,repro-openmp -r rts --overwrite --no-transfer -s -e SPEAR_c192_o1_Hist_AllForc_IC1921_Q50_ens_01_03
multi-frerun -x SPEAR.xml -p ncrc5.intel23-classic -t repro-openmp,prod-openmp -r rts --overwrite --no-transfer -s -e SPEAR_c192_o1_Control_1850_Q50
```

```
frecheck -x SPEAR.xml -p ncrc5.intel23-classic -t prod-openmp -r rts SPEAR_c192_o1_Control_1850_Q50
frecheck -x SPEAR.xml -p ncrc5.intel23-classic -t repro-openmp -r rts SPEAR_c192_o1_Control_1850_Q50
frecheck -x SPEAR.xml -p ncrc5.intel23-classic -t repro-openmp -r rts SPEAR_c192_o1_Hist_AllForc_IC1921_Q50_ens_01_03
frecheck -x SPEAR.xml -p ncrc5.intel23-classic -t prod-openmp -r rts SPEAR_c192_o1_Hist_AllForc_IC1921_Q50_ens_01_03
```

Canopy Bare Metal Testing:
```
module load fre/2025.01
fre make run-fremake -y spear.yaml -p ncrc5.intel23-classic -t prod-openmp --execute

module load fre/bronx-23
frerun -x SPEAR.xml -p ncrc5.intel23-classic -t prod-openmp -r rts --overwrite --no-transfer -s SPEAR_c192_o1_Control_1850_Q50_canopy
frerun -x SPEAR.xml -p ncrc5.intel23-classic -t prod-openmp -r rts --overwrite --no-transfer -s SPEAR_c192_o1_Hist_AllForc_IC1921_Q50_ens_01_03_canopy

frecheck -x SPEAR.xml -p ncrc5.intel23-classic -t prod-openmp -r rts SPEAR_c192_o1_Control_1850_Q50_canopy
frecheck -x SPEAR.xml -p ncrc5.intel23-classic -t prod-openmp -r rts SPEAR_c192_o1_Hist_AllForc_IC1921_Q50_ens_01_03_canopy
```

Canopy Container Testing:
```
module load fre/2025.01
setenv TMPDIR /tmp/containers/$USER
fre make run-fremake -y spear.yaml -p hpcme.2023 -t prod-openmp -npc --execute

module load fre/bronx-23
frerun -x SPEAR.xml -p ncrc5.intel23-classic -t prod-openmp -r rts --unique --no-transfer -s SPEAR_c192_o1_Control_1850_Q50_canopy --container
frerun -x SPEAR.xml -p ncrc5.intel23-classic -t prod-openmp -r rts --unique --no-transfer -s SPEAR_c192_o1_Hist_AllForc_IC1921_Q50_ens_01_03_canopy --container

```

