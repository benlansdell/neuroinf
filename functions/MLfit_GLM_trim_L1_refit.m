function [gg, fval,H] = MLfit_GLM_trim_refit(gg,Stim,optimArgs,processed, trim, offset);
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

thresh = 1e-5;
% Set optimization parameters 
if nargin > 2
    opts = optimset('Gradobj','on','Hessian','on', optimArgs{:});
else
    opts = optimset('Gradobj','on','Hessian','on','display','iter');
end
%Check which coupling terms are zero and mask these ones out, then redo the fitting with just MLE
%Anything less than 10^-5 counts as zero
gg.mask = any(abs(gg.ih2) > thresh,1);
% Set initial params
prs1 = extractFitPrs_GLM_trim(gg,Stim,MAXSIZE,processed, trim, offset);
% minimize negative log likelihood. The warm-up solution for L1 reg
[prs,fval] = fminunc(@(p) Loss_GLM_logli(p, gg.mask), prs1, opts);
if nargout > 2 % Compute Hessian if desired
    [fv,gradval,H] = Loss_GLM_logli(prs);
end
% Put returned vals back into param structure ------
gg = reinsertFitPrs_GLM(gg,prs);