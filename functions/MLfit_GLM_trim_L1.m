function [gg, fval,H] = MLfit_GLM_trim(gg,Stim,optimArgs,processed, trim, offset, lambda);
%  [ggnew,fval,H] = MLfit_GLM(gg,Stim,optimArgs);
% 
%  Computes the ML estimate for GLM params, using grad and hessians.
%  Assumes basis for temporal dimensions of stim filter
%
%  Inputs: 
%     gg = param struct
%     Stim = stimulus
%     optimArgs = cell array of optimization params (optional)
%
%  Outputs:
%     ggnew = new param struct (with ML params);
%     fval = negative log-likelihood at ML estimate

%Only include times within stimulus that are within trial

MAXSIZE  = 1e7;  % Maximum amount to be held in memory at once;
if (nargin < 6) offset = 1; end
if (nargin < 7) lambda = 1; end

% Set optimization parameters 
if nargin > 2
    opts = optimset('Gradobj','on','Hessian','on', optimArgs{:});
else
    opts = optimset('Gradobj','on','Hessian','on','display','iter');
end

% Set initial params
prs0 = extractFitPrs_GLM_trim(gg,Stim,MAXSIZE,processed, trim, offset);

% minimize negative log likelihood. The warm-up solution for L1 reg
[prs,fval] = fminunc(@Loss_GLM_logli,prs0,opts);
if nargout > 2 % Compute Hessian if desired
    [fv,gradval,H] = Loss_GLM_logli(prs);
end

groups = zeros(size(prs));
nG = size(gg.tsp2,2);
nP = size(gg.ihbas2,2);
gOptions.maxIter = 2000;
gOptions.verbose = 2; % Set to 0 to turn off output
gOptions.corrections = 10; % Number of corrections to store for L-BFGS methods
gOptions.norm = 2; % Set to inf to use infinity norm
options = gOptions;
istart = size(gg.kt,1)*size(gg.kt,2)+1+size(gg.ih,1);
for idx = 1:nG
	indices = istart + (((idx-1)*nP+1):(idx*nP));
	groups(indices) = idx; 
end
lambdaVect = lambda*ones(nG, 1);	
w = L1GeneralGroup_Auxiliary(@Loss_GLM_logli,prs,lambdaVect,groups,options);

% Put returned vals back into param structure ------
gg = reinsertFitPrs_GLM(gg,w);

%----------------------------------------------------
% % ------ Check analytic gradients, Hessians -------
% HessCheck(@Loss_GLM_logli,prs0,opts);
% HessCheck_Elts(@Loss_GLM_logli, [1 12],prs0,opts);
% tic; [lival,J,H]=Loss_GLM_logli(prs0); toc;