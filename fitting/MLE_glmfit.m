function model = MLE_glmfit(data, const)
	%Fit GLM using MATLAB's glmfit
	%     
	%Input:
	%	data = covariate data output structure from any function in ./models
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
	%		converged = [nU x 1] array listing 1 if the IRLS converged within iteration limit, 0 if not
	%		conditioned = [nU x 1] array listing 1 if the IRLS did not issue an ill-conditioned warning, 0 if it did
	%
	%Test code:
	%	responsefiles = dir('./data/*.isk');
	%	rf = responsefiles(1:3);
	%	stimfile = './data/whitenoise.raw';
	%	binsize = 1;
	%	unit = 10;
	%	nK_sp = 20;
	%	nK_stm = 6;
	%	const = 'on';
	%	processed = preprocess(stimfile, rf, binsize, unit);
	%	data = filters_sprc_stm(processed, nK_sp, nK_stm);
	%	model = MLE_glmfit(data, const);

	if (nargin < 2) const = 'on'; end

	nU = 1;
	nK = size(data.X,2);
	if strcmp(const, 'on')
		model.b_hat = zeros(nU, nK+1);
	else
		model.b_hat = zeros(nU, nK);
	end
	model.dev = cell(nU,1);
	model.stats = cell(nU,1);
	model.converged = ones(nU,1);
	model.conditioned = ones(nU,1);
	%For each unit, fit a GLM to the torque data
	display(['Fitting GLM by MLE with IRLS. Fitting ' num2str(nU) ' units.'])
	[b, dev, stats] = glmfit_matlab(squeeze(data.X),data.y(1,:),'poisson', 'constant', const);
	%Catch if a warning was raised about badly conditioned matrix
	[warn, warnid] = lastwarn;
	if ~strcmp(warn, '')
   		switch warnid
       	case 'stats:glmfit:IterationLimit'
       		model.converged(1) = 0;
       	case 'stats:glmfit:BadScaling'
       		model.conditioned(1) = 0;
      		end
    end
    lastwarn('')
	%Extract filters fitted...
	model.b_hat(1,:) = b;	
	model.dev{1} = dev;
	%Remove residual components since these take up a lot of memory
	model.N = size(stats.resid,1);
	stats = rmfield(stats, {'resid', 'residp', 'residd', 'resida', 'wts'});
	model.stats{1} = stats;
	model.nspikes(1) = sum(data.y(1,:));

	%model.logli = ll(model, data, 'poisson');
	if ~strcmp(const, 'on')
		model.b_hat = [zeros(nU, 1), model.b_hat]
	end
	display('Done')