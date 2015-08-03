function model = L1_grouplasso(processed, data, lambda, const)
	%Fit GLM using fminunc
	%	Model includes a penalty term for complexity in the network filters using
	%	group LASSO
	%     
	%Input:
	%	data = covariate data output structure from any function in ./models
	%	lambda = weight for filter complexity penalty
	%	const = (optional, default = 'on') whether to fit a constant term to the model or not, 'on' or 'off'
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
	%	datafile = './data/mabel_reaching_5-4-10.mat';
	%	binsize = 1/100;
	%	nK_sp = 20;
	%	nK_stm = 6;
	%	dt_sp = binsize;
	%	dt_stm = 5/100;
	%	lambda = 10;
	%	const = 'on';
	%	unitidx = 13;
	%	processed = preprocess_monkey(datafile, binsize, unitidx);
	%	data = filters_monkey_sp_stm_network(processed, nK_sp, nK_stm, dt_sp, dt_stm);
	%	model_sp = L1_grouplasso(processed, data, lambda, const);

	if (nargin < 3) const = 'on'; end

	%Truncate for testing...
	data.X = data.X(1:10000,:);
	data.y = data.y(1:10000);
	%Make smaller groups for testing...
	data.k = data.k(1:13,:);
	data.X = [data.X(:,1:260)];

	nU = size(data.y,1);
	nK = size(data.X,2);
	nG = size(data.k,1);
	nB = size(data.X,1);

	%Warm up solution with MLE
	model = MLE_SD_grouplasso(data, lambda, const);
	unitidx = processed.unitidx;
	if strcmp(const, 'on')
		groups = zeros(nK+1,1);
		%model.b_hat = zeros(nU, nK+1);
		data.X = [ones(nB, 1), data.X];
	else
		%model.b_hat = zeros(nU, nK);
		groups = zeros(nK,1);
	end
	model.dev = cell(nU,1);
	model.stats = cell(nU,1);

	b0 = model.b_hat;
	%Need a better intial guess if this is to run well...
	%b0 = 1e-5*ones(size(model.b_hat))';
	gOptions.maxIter = 2000;
	gOptions.verbose = 2; % Set to 0 to turn off output
	gOptions.corrections = 10; % Number of corrections to store for L-BFGS methods
	gOptions.norm = 2; % Set to inf to use infinity norm
	options = gOptions;
	funObj = @(W) logli(data.X, data.y, W);
	for idx = 1:nG
		indices = data.k{idx,2};
		if strcmp(const, 'on')
			indices = indices + 1;
		end
		%Don't apply sparsity to own spike history filter
		if idx ~= unitidx
			groups(indices) = idx; 
		else
			groups(indices) = 0;
		end
	end
	lambda = 1e-2;
	lambdaVect = lambda*ones(nG-1, 1);	
	w = L1GeneralGroup_Auxiliary(funObj,b0',lambdaVect,groups,options);

	options.method = 'opg';
	options.L = 1/1000;
	w = L1GeneralGroup_Auxiliary(funObj,b0',lambdaVect,groups,options);

	fprintf('\nProjected Quasi-Newton\n');
	options = gOptions;
	options.method = 'pqn';
	w = L1GeneralGroup_Auxiliary(funObj,b0',lambdaVect,groups,options);
	Wpqn = reshape(w,nVars,nTargets);
	pause;

	fprintf('\nBarzilai-Borwein Soft-Threshold\n');
	options = gOptions;
	w = L1GeneralGroup_SoftThresh(funObj,b0',lambdaVect,groups,options);
	Wbbst = reshape(w,nVars,nTargets);
	pause;

	fprintf('\nQuasi-Newton Soft-Threshold\n');
	options = gOptions;
	options.method = 'qnst';
	w = L1GeneralGroup_SoftThresh(funObj,b0',lambdaVect,groups,options);
	Wbbst = reshape(w,nVars,nTargets);

	lambdas = [1e-1 3e-1 1];
	for lambda = lambdas
		lambdaVect = lambda*ones(nG-1, 1);
		w = L1GeneralGroup_Auxiliary(funObj,w,lambdaVect,groups,options);
	end		

	model.b_hat = w;
	if ~strcmp(const, 'on')
		model.b_hat = [zeros(1, 1), model.b_hat];
	end

	model.nspikes(1) = sum(data.y(1,:));
end

function [ll, d, h] = logli(X,y,b_hat)
	ll = -l(X,y,b_hat);
	d = -dl(X,y,b_hat);
	h = -hessian(X,y,b_hat);
end

function ll = l(X,y,b_hat,r,lambda)
	ll = y*X*b_hat-sum(exp(X*b_hat)+log(y+0.00001)');
end

function d = dl(X,y,b_hat)
	d = (y*X)'-X'*exp(X*b_hat);
end

function h = hessian(X,y,b_hat)
	nT = size(X,1);
	h = -X'*spdiags(exp(X*b_hat), 0, nT, nT)*X;
end