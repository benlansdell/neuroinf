function model = MLE_SD(data, const)
	%Fit GLM to spike data using MATLAB's fminunc
	%     
	%Input:
	%	data = covariate data output structure from any function in ./models
	%	const = (optional, default = 'on') whether to fit a constant term to the model 
	%		or not, 'on' or 'off'
	%   
	%Output:
	%	model is a structure containing the following fields:
	%		b_hat = [nU x (nK + 1)] array of estimated filter coefficients. Refer
	%			to input 'data' structure for identity of each param.
	%			Note: if a constant term is not fit, a column of zeros is appended 
	%			to b_hat to make dimensions consistent
	%		dev = [nU x 1] cell array listing deviance of each unit's fit
	%		stats = [nU x 1] cell array listing fitting statistics output from glmfit
	%
	%Test code:
	%	responsefiles = dir('./data/*.isk');
	%	rf = responsefiles(1:3);
	%	stimfile = './data/whitenoise.raw';
	%	binsize = 1;
	%	unit = 10;
	%	nK_sp = 6;
	%	nK_stm = 6;
	%	const = 'on';
	%	processed = preprocess(stimfile, rf, binsize, unit);
	%	data = filters_sp_stm_network(processed, nK_sp, nK_stm);
	%	model_sp = MLE_SD(data, const);

	if (nargin < 2) const = 'on'; end
	model = MLE_SD_reg(data, 0, const);
end