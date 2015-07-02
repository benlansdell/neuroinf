function [model, intermediates] = MLE_SD_rc_reg(data, lambda, const)
	%Fit GLM to spike data from blackrock recording file for each unit above a specified threshold using fminunc
	%	Model includes a penalty term for complexity in the network filters
	%	Weighted by parameter lambda.
	%	Assume stim filters use a raised cosine basis
	%     
	%Input:
	%	data = covariate data output structure from any function in ./models
	%	lambda = weight for 
	%	const = (optional, default = 'on') whether to fit a constant term to the model or not, 'on' or 'off'
	%   
	%Output:
	%	model is a structure containing the following fields:
	%		b_hat = [nU x (nK + 1)] array with spikes from all channels binned according to binsize. nB = no. bins, nU = no. units.
	%			Note: if a constant term is not fit, a column of zeros is appended to b_hat to make dimensions consistent
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
	%	lambda = 0.1;
	%	processed = preprocess(stimfile, rf, binsize, unit);
	%	data = filters_sprc_stm_network(processed, nK_sp, nK_stm);
	%	model_sp = MLE_SD_rc_reg(data, lambda, const);

	if (nargin < 3) const = 'on'; end

	nU = size(data.y,1);
	nK = size(data.X,2);
	if strcmp(const, 'on')
		model.b_hat = zeros(nU, nK+1);
	else
		model.b_hat = zeros(nU, nK);
	end
	model.dev = cell(nU,1);
	model.stats = cell(nU,1);

	%Need a better intial guess if this is to run well...
	b0 = 1e-5*ones(size(model.b_hat))';
	%b0 = zeros(size(model.b_hat))';

	%Regularization parameters
	nC = size(data.spbasis,1);
	nRC = size(data.rcbasis,1);
	%Assume all but the last model component is a spike history term...
	cplidx = [];
	for idx = 1:(size(data.k,1)-1)
		cplidx = [cplidx(:); data.k{idx,2}(:)];
	end

	%Fit GLm
	[b, dev] = glmfit_SD_reg(data.X,data.y(1,:), const, b0, nC, nRC, cplidx, data.spbasis, lambda);
	%Extract filters fitted...
	model.b_hat(1,:) = b;	
	model.dev{1} = dev;
	model.stats{1} = 0;
	model.nspikes(1) = sum(data.y(1,:));

	if ~strcmp(const, 'on')
		model.b_hat = [zeros(1, 1), model.b_hat];
	end

	%TODO: report convergence
	model.conditioned = 1;
	model.converged = 1;
end

function [b, dev] = glmfit_SD_reg(X, y, const, b0, nC, nRC, cplidx, spbasis, lambda)
	%X is a matrix of data for that unit
	%y is the set of observations for that unit
	if strcmp(const, 'on')
		X = [ones(size(X,1),1), X];		
	end

	db = 0.1;
	N = size(X, 1);
	nIter = 5e7;
	nB = size(X, 2);
	tolfun = 1e-9;
	opts = optimset('Display', 'iter', 'MaxIter', 4000, 'GradObj', 'on', 'Hessian', 'on', 'TolFun', tolfun);
	[b, fval, exitflag, output] = fminunc(@(b) logli(X,y,b,cplidx,nC,nRC,spbasis,lambda), b0, opts);

	%TODO compute deviance
	dev = 0;
end

function [ll, d, h] = logli(X,y,b_hat,cplidx,nC,nRC,spbasis,lambda)
	r = reshape(b_hat(cplidx+1),nRC,[]);
	ll = -l(X,y,b_hat,r,spbasis,lambda);
	d = -dl(X,y,b_hat,r,spbasis,cplidx,lambda);
	h = -hessian(X,y,b_hat,r,spbasis,cplidx,lambda);
end

function ll = l(X,y,b_hat,r,spbasis,lambda)
	%Coupling filters in standard basis
	cpl = spbasis*r;
	%Log likelihood
	ll = y*X*b_hat-sum(exp(X*b_hat)+log(y+0.00001)');
	%Add regularization term (we're taking the one norm of a two norm...)
	ll = ll-lambda*sum(sqrt(sum(cpl.^2,2)),1);
	%Note that doing this
	%ll = ll- lambda*sum(sqrt(sum(cpl.^2,1)),2);
	%instead gives a group LASSO-like penalty
end

function d = dl(X,y,b_hat,r,spbasis,cplidx,lambda)
	%Coupling filters in standard basis
	cpl = spbasis*r;
	%Derivative of log likelihood with respect to each covariate b
	d = (y*X)' - X'*exp(X*b_hat);
	%Regularization term
	dp = -lambda*spbasis'*diag(sqrt(sum(cpl.^2,2)))*cpl;
	%Reshape matrix of filters back into vector of parameters
	dp = reshape(dp,1,[])';
	d(cplidx+1) = d(cplidx+1)+dp;
end

function h = hessian(X,y,b_hat,r,spbasis,cplidx,lambda)
	nB = length(cplidx);
	nRC = size(r,1);
	nF = size(r,2);
	nT = size(X,1);
	h = -X'*spdiags(exp(X*b_hat), 0, nT, nT)*X;
	%Coupling filters in standard basis
	cpl = spbasis*r;
	%The less efficient way
	%nB = length(b_hat);
	%h = zeros(nB);
	%for i = 1:nB
	%	for j = 1:nB
	%		h(i,j) = -sum(X(:,i).*X(:,j).*exp(X*b_hat));
	%	end
	%end
	%Regularization, much slower...
	hp = zeros(nRC, nF, nRC, nF);
	for kp = 1:nRC
		for jp = 1:nF
			for lp = 1:nRC
				for mp = 1:nF
					hp(kp,jp,lp,mp) = -lambda*(-sum(spbasis(:,kp).*spbasis(:,lp).*cpl(:,jp).*cpl(:,mp).*(power(sum(cpl.^2,2),-1.5))));
					if jp == mp
						hp(kp,jp,lp,mp) = hp(kp,jp,lp,mp)-lambda*(sum(spbasis(:,kp).*spbasis(:,lp).*power(sum(cpl.^2,2),-0.5)));
					end
				end
			end
		end
	end
	%Reshape and add to hessian matrix
	hp = reshape(hp,nRC,nF,1,[]);
	hp = reshape(hp,1,[],1,nB);
	hp = squeeze(hp);
	h((cplidx+1),(cplidx+1)) = h((cplidx+1),(cplidx+1))+hp;
end