lambdas = [.1 .3 1 3 10 30 100 300];
wd = '../MotorData/';

%Simulate coupled GLM with range of lambda values
%This takes a long time, so split into many runs of a smaller number of
%repetitions to divide between CPUs/computers if necessary. 

for idx = 1:length(lambdas)
	sim_coupled_GLM_L1(wd, 1, 2, idx);
	sim_coupled_GLM_L1(wd, 2, 2, idx);
	sim_coupled_GLM_L1(wd, 3, 2, idx);
	sim_coupled_GLM_L1(wd, 4, 2, idx);
	sim_coupled_GLM_L1(wd, 5, 2, idx);
	sim_coupled_GLM_L1(wd, 6, 2, idx);
	sim_coupled_GLM_L1(wd, 7, 2, idx);
	sim_coupled_GLM_L1(wd, 8, 2, idx);
	sim_coupled_GLM_L1(wd, 9, 2, idx);
end