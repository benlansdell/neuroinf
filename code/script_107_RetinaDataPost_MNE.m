% This script computes the Maximum Noise Entropy (MNE) LN model using a 
% gradient ascent algorithm. In contrast to the STA and STC models, here 
% the nonlinearity and the relevant subspace are computed together. 

% Fitting is a convex problem with respect to the the training set, meaning
% that there are no local minima. However, convergence on the test set 
% often leads to over-fitting. Therefore the trainig set is split into 4 
% parts and fitting is repeated 4 times: 
% in each repeated fitting ('jackknife') the gradient ascent runs and the 
% likelihood function is computed on the training (3/4) and test (1/4) 
% parts of the data. Because the problem is convex, the likelihood of the
% training portion always increases. The likelihood function computed on 
% the test portion starts to decrease when the model is over-fit, so
% the algorithm stops at that point.
% Included Null hypothesis calcualtion.


% Averaging and significance of the J matrix - plus generates predictions


clear ; 
cd('/Users/Shared/MANUSCRIPTS/Adrienne_Yonitan/Retinal_data') %DK files
%cd('/home/aljadeff/Documents/MATLAB/neuroinfo_retina/data') %YA files

stim_length = {'short','long'} ;

sr = 30 ;          % (Hz) sampling rate
dt = 1/sr ;        % delta t of stimulus 

nspk = zeros(1,53) ;
logl_mne = zeros(2,53) ;

order   = 2 ;   % order of MNE model to fit
njack   = 4 ;   % # jackknives to run (also determines the size of each jackknives)

rep = 500 ;        % number of repeats used to compute the null distribution of covariance eigenvalues

a = 0 ;

for icell = 1:53
    for iL = 1:2
        disp(num2str(icell)) ; disp(num2str(iL)) ; % counter as loop is long

        load(['Retina_cell_' num2str(icell) '_stim_resp_' stim_length{iL} '.mat']) ;
        load(['Retina_cell_' num2str(icell) '_mne_model_' stim_length{iL} '.mat']) ; 
    
        S = S - repmat(mean(S),[T,1]);
        S = S./repmat(std(S),[T,1]);
        R = sign(R) ;
    
        St = St - repmat(mean(St),[Tt,1]);
        St = St./repmat(std(St),[Tt,1]);
        Rt = sign(Rt) ;
    
        m = mean(mne_model,1) ;
        A = m(1) ;
        H = m(2:p+1) ;
        J = reshape(m(p+2:end),p,p) ;
    
        [vJ,eJ] = eig(J) ; 
        [eJ,iJ] = sort(diag(eJ)) ;
        vJ = vJ(:,iJ) ;
    
        eJnull = zeros(1,p*rep) ;     
    
        dJ = diag(J) ;
        oJ = reshape(triu(J,1),1,p*p) ;
        [~,iJ] = sort(abs(oJ),'descend') ;
        iJ = iJ(1:p*(p-1)/2) ;
        oJ = oJ(iJ) ;
    
        for ir = 1:rep
            dinull = randperm(p) ;
            oinull = randperm(p*(p-1)/2) ;
            Jnull  = zeros(p) ;
            Jnull(iJ(oinull)) = oJ ;
            Jnull = Jnull + Jnull' + diag(dJ(dinull)) ;
            eJnull((ir-1)*p+(1:p)) = eig(Jnull) ;
        end
    
        max_enull = prctile(eJnull(:),100-a) ;  % upper bound of null distribution  
        min_enull = prctile(eJnull(:),a)  ;     % lower bound of null distribution
    
        isig = [find(eJ>max_enull) ; find(eJ<min_enull)] ;
        nsig = length(isig) ;                   % number of significant quadratic features
    
         if nsig>0 
            [~,iis] = sort(abs(eJ(isig)),'descend') ;
            isig = isig(iis) ;                   % orders the quadratic features according to the absolute value of their corresponding eigenvalues   
         end
   
        save(['Retina_cell_' num2str(icell) '_mne_' stim_length{iL} '.mat'],'mne_model','p','pT','pX','nsig','J','H','A','isig','eJ','vJ','eJnull','max_enull','min_enull','rep') ;
    
    SHt  = St*H' ;
    SJSt = sum(St.*(St*J),2) ;

    Rt_mne = 1./(1+exp(SJSt+SHt+A)) ;
    logl_mne(iL,icell) = mean(Rt.*log(Rt_mne)-Rt_mne*dt) ;
    end
end

save('Retina_MNE_model_All_LogLikelihood.mat','logl_mne') ;