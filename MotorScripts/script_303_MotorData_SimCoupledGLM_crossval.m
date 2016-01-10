lambdas = [.1 .3 1 3 10 30 100 300];
wd = '../MotorData/';

%Simulate coupled GLM with range of lambda values
%This takes a long time, so split into many runs of a smaller number of
%repetitions to divide between CPUs/computers if necessary. 

%These runs should be moved to the one directory, and combined and analyzed afterwards
%by collapse_sim_coupled_GLM_L1.m


nfolds = 5;
running = [1];
for fold = running
	for idx = 1:length(lambdas)
		sim_coupled_GLM_L1_crossval(wd, 1, 83, idx, fold, nfolds);
		sim_coupled_GLM_L1_crossval(wd, 2, 83, idx, fold, nfolds);
		sim_coupled_GLM_L1_crossval(wd, 3, 83, idx, fold, nfolds);
	end
end

for fold = running
	for idx = 1:length(lambdas)
		sim_coupled_GLM_L1_crossval(wd, 4, 83, idx, fold, nfolds);
		sim_coupled_GLM_L1_crossval(wd, 5, 83, idx, fold, nfolds);
		sim_coupled_GLM_L1_crossval(wd, 6, 83, idx, fold, nfolds);
	end
end

for fold = running
	for idx = 1:length(lambdas)
		sim_coupled_GLM_L1_crossval(wd, 7, 83, idx, fold, nfolds);
		sim_coupled_GLM_L1_crossval(wd, 8, 83, idx, fold, nfolds);
		sim_coupled_GLM_L1_crossval(wd, 9, 83, idx, fold, nfolds);
	end
end
%Cross validate simulations by only taking subsamples of data...