
% This script computes the Spike Triggered Average (STA) and the STA LN
%  model (i.e. the 1 dimensional linear nonlinear model where the single
%  dimension - or feature - is the STA). The nonlinearity is computed using
%  Bayes' rule as explained in the manuscript. 

% The STA model is computed in three steps: 
% 1. computation of the spike triggered avereage (i.e. reducing the 
%    dimensionality of the neuron's stimulus dependence to 1). 
% 2. computation of the spiking nonlinearity (using Bayes' rule and binned 
%    probability distributions)
% 3. Building an LN for general stimulus by interpolating the spiking 
%    nonlinearity.
    

% Given the 1d LN model and the test data set, the log likelihood of the
%  model is also computed. 

clear ; 
%cd('/Users/Shared/MANUSCRIPTS/Adrienne_Yonitan/Retinal_data') %DK files
%cd('/home/aljadeff/Documents/MATLAB/neuroinfo_retina/data') %YA files
cd('/home/lansdell/projects/neuroinf/data') %BL files

stim_length = {'short','long'} ;

dt = 1/30 ;%ds/sr ;    % delta t of stimulus 

nspk = zeros(1,53) ;       % a vector with the number of spikes for each cell 
logl_sta   = zeros(53,2) ; % a vector with the log likelihood of the STA model for each cell 
logl_0   = zeros(53,2) ;   % a vector with the null value of the log likelihood for each cell, 
                          % assuming each cell fires at random with a probability fixed by the average firing rate found in the experiment

nbns = 20 ;               % number of bins used to construct the prior, conditional and posterior probability distributions
%tp = -dt*(pT:-1:1) ;       % the time vector associated with the stimulus history

linearrect = @(x) x.*(sign(x)+1)/2 ;

for icell = 1:53
    for iL = 1:2
        load(['Retina_cell_' num2str(icell) '_stim_resp_' stim_length{iL} '.mat']) ;
        iR = find(R) ; 
        
        nspk(icell) = sum(R) ;

        % 1. computation of the spike triggered avereage:
        % -----------------------------------------------
    
        sta = S'*R/sum(R) - mean(S,1)' ;    % the 1st term on the RHS is the stimulus average conditioned on a spike
                                            % the 2nd term on the RHS is the unconditioned stimulus average 
                                            % here the STA is a (1 x p) vector
    
        sta = sta / norm(sta) ;             % the STA is transposed to give a (p x 1) vector and normalized to 1 
    
    
        % 2. computation of the spiking nonlinearity 
        % ------------------------------------------
    
        Ps = nspk(icell)/T ;              % Ps  - probability of a spike                         
                                          %     - the number of spikes divided by the number of stimulus 'frames'    
        [Pf,ctrs] = hist(S*sta,nbns) ;    % Pf  - probability of a stimulus projected on feature 
                                          %     - a normalized count of the projections of stimuli on the STA in nbns bins  
        dbns = ctrs(2)-ctrs(1) ;          % the bin size
        Pf   = Pf/T/dbns ;                % normalizing Pf to 1   
        Pfs = hist(S(iR,:)*sta,ctrs) ;    % Pfs - probability of a stimulus projected on feature given a spike
                                          %     - a normalized count of the projections of stimuli on the STA in nbns bins conditioned on a spike 
        Pfs = Pfs/nspk(icell)/dbns ;          % normalizing Pfs to 1
        Psf    = Pfs./Pf*Ps ;             % Psf - probability of a spike given stimulus projected on a feature
                                          %     - computed by invoking Bayes' rule (see details in text)
        Psf(isnan(Psf)) = 0 ;             % setting NaNs to 0 
    
        Psf = Psf' ;                      % transposing Psf to give a (nbns x 1) vector
        
        % 3. Building an LN for general stimulus
        % --------------------------------------
        
        sta_model = @(x)interp1(ctrs,Psf,x,'pchip') ;                       % first the STA model is defined to be a 
                                                                            % smoothing spline interpolation of the spiking
                                                                            % nonlinearity computed using binned probability distributions                                
        sta_model_rect = @(x) linearrect(sta_model(x))+1e-8 ;  % then the model passed through a linear rectifier, to 
                                                                            % assure that the probability of a spike is non-negative
                                                                            % the infinitesimal constant offset (1e-8) ensures that the 
                                                                            % loglikelihood calculation below does not diverge
    
        Ps_model = sum(sta_model_rect(S*sta)) ;                             % finally the model has to be renormalized such that the number of 
                                                                            % of predicted spikes (on the training set) is fixed to what was 
                                                                            % found in the experiment. This step is needed because of the 
                                                                            % smoothing and rectifying operations  
        sta_model_rect_norm = @(x) sta_model_rect(x)*sum(iR)/Ps_model ;     % this normalization is done by computing the number of predicted 
                                                                            % spikes from the smooth rectified model (Ps_model)                                                                 %                                                                         
                                                                            % then the model multiplied by the ratio of the actual number of spikes and Ps_model                                                                        % 
    
        % Log likelihood                                                                    
        logl_sta(icell,iL) = mean(Rt.*log(sta_model_rect_norm(St*sta))-sta_model_rect_norm(St*sta)/dt) ;
        logl_0(icell,iL)   = mean(Rt.*log(sum(Rt/Tt))-Rt/(Tt*dt)) ;
        save(['Retina_cell_' num2str(icell) '_sta_' stim_length{iL} '.mat'],'sta_model_rect_norm','sta','pT','pX','p','ctrs','Pfs','Psf','Pf') ; 
    end
end
save('Retina_STA_model_All_LogLikelihood.mat','logl_sta','logl_0','nspk') ;