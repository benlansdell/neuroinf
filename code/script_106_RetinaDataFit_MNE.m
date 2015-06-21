% This script computes the Maximum Noise Entropy (MNE) LN model using a 
% gradient ascent algorithm. In contrast to the STA and STC models, here 
% the nonlinearity and the relevant subspace are computed together. 

% Fitting is a convex problem with respect to the the training set, meaning
% that there are no local minima. However, convergence on the test set 
% often leads to over-fitting. Therefore the trainig set is split into 4 
% parts and fitting is repeated 4 times.
% In each repeated fitting ('jackknife') the gradient ascent runs and the
% likelihood function is computed on the training (3/4) and test (1/4) 
% parts of the data. Because the problem is convex, the likelihood of the
% training portion always increases. The likelihood function computed on 
% the test portion starts to decrease when the model is over-fit, so
% the algorithm stops at that point.

% Fit's the J matrix


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

for icell = 1:53
    for iL = 1:2
        disp(num2str(icell)) ; disp(num2str(iL)) ; % counter as loop is long

        load(['Retina_cell_' num2str(icell) '_stim_resp_' stim_length{iL} '.mat']) ;
    
        S = S - repmat(mean(S),[T,1]);
        S = S./repmat(std(S),[T,1]);
        R = sign(R) ;  
        
        St = St - repmat(mean(St),[Tt,1]);
        St = St./repmat(std(St),[Tt,1]);
        Rt = sign(Rt) ;
        
        mne_model = zeros(njack,p*p+p+1) ;

        for j = 1:njack
            Tj = floor(T/njack);  % could be changed to do different # of jackknives
            jtst = (1+(j-1)*Tj:j*Tj) ; 
            jfit = setdiff(1:T,jtst) ;
            Sjtst = S(jtst,:) ;
            Rjtst = R(jtst) ;
            Sjfit  = S(jfit,:) ;
            Rjfit  = R(jfit) ;
            mne_model(j,:) = MNEfit_RetinaData(Sjfit, Rjfit, Sjtst, Rjtst, order) ;
        end
        save(['Retina_cell_' num2str(icell) '_mne_model_' stim_length{iL} '.mat'],'mne_model','p') ;
    end
end