#!/bin/sh
cp ./scripts/* ./mcc
cp ./functions/* ./mcc
cp ./GLM_Algorithm_Functions/* ./mcc
cp ./GLM_Algorithm_Functions/tools_mexcode/* ./mcc
cp ./GLM_Algorithm_Functions/tools_misc/* ./mcc
cp ./GLM_Algorithm_Functions/tools_splines/* ./mcc
cp ./GLM_Algorithm_Functions/nlfuns/* ./mcc
cd mcc
matlab -nojvm -nodisplay -nosplash -r "compilesimcode"
cp sim_coupled_GLM bin
cp sim_coupled_GLM_L1 bin
cp sim_coupled_GLM_L1_MLE bin
cp run_sim_coupled_GLM.sh bin
cp run_sim_coupled_GLM_L1.sh bin
cp run_sim_coupled_GLM_L1_MLE.sh bin
cd ..
