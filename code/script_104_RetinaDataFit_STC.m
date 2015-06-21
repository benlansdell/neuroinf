% This script computes the Spike Triggered Covariance (STC) and the STC LN
%  model (i.e. the multi-dimensional linear nonlinear model where the dimesnions 
%  - or features - are the STA and STC). The nonlinearity is computed using
%  Bayes' rule as explained in the manuscript. 

% The STC model is computed in three steps: 
% 1. computation of the spike triggered average, and projection of STA out of stimulus.
% 2. computation of the spike triggered covariance (Cs) and the matrix of covariance differences (dC) 
% 3. finding significant STC dimensions by comparing to a null eigenvalue distribution 
% 4. computation of the spiking nonlinearity (using Bayes' rule and binned probability distributions)
% 5. building an LN for general stimulus by interpolating the spiking nonlinearity.
    
% Given the LN model and the test data set, the log likelihood of the model is also computed. 


clear ; 
cd('/Users/Shared/MANUSCRIPTS/Adrienne_Yonitan/Retinal_data') %DK files
%cd('/home/aljadeff/Documents/MATLAB/neuroinfo_retina/data') %YA files

stim_length = {'short','long'} ;

dt = 1/30 ;%ds/sr ;    % delta t of stimulus 

nspk = zeros(1,53) ;       % a vector with the number of spikes for each cell 

rep = 1000 ;        % number of repeats used to compute the null distribution of covariance eigenvalues

a = 0 ;            % level of significance required to call an STC feature significant
                   % if an eigenvalue of dC is smaller than the a-th percentile of the null eigenvalue 
                   % distribution or larger than the 100-a percentile, it will be determined to be significant

for icell = 1:53 ;
    for iL = 1:2 ;
      disp(num2str(icell)) ; disp(num2str(iL)) ; % counter as loop is long

      load(['Retina_cell_' num2str(icell) '_stim_resp_' stim_length{iL} '.mat']) ;
        iR = find(R) ; 
        
        nspk(icell) = sum(R) ;

        % 1. computation of the spike triggered avereage:
        % -----------------------------------------------
    
        sta = S'*R/sum(R) - mean(S,1)' ;    % the 1st term on the RHS is the stimulus average conditioned on a spike
                                            % the 2nd term on the RHS is the unconditioned stimulus average 
                                            % here the STA is a (1 x p) vector
    
        sta = sta / norm(sta) ;             % the STA is transposed to give a (p x 1) vector and normalized to 1 
        
        S1 = S-(S*sta)*sta' ;               % S1 is the stimulus with the STA projected out
                                            % meaning that every stimulus vector in S1 is 
                                            % orthogonal to the STA but has the same projections 
                                            % on the rest of the stimulus dimensions

        Cp = cov(S1-repmat(mean(S1,1),[T,1])) ;  % Cp is the prior covariance matrix, 
                                             % meaning the covariance of mean subtracted and stimulus 
                                             % with the STA projected out
                                             
        Cs = 1/mean(R)*cov(S1.*repmat(sqrt(R),[1 p])-repmat(mean(S1.*repmat(sqrt(R),[1 p])),[T 1])) ; 
                                             % Cs is the spike triggered covariance matrix, 
                                             % meaning the covariance of mean subtracted and stimulus 
                                             % with the STA projected out, conditioned on a spike
                                             
        dC = Cs-Cp ;                         % the STC dimensions are found by diagonalizing 
        [vdC,edC] = eig(dC) ;                % dC - the matrix of covariance differences    
        [edC,idC] = sort(diag(edC)) ;
        vdC = vdC(:,idC) ;
    
        edCnull = zeros(1,p*rep) ;           % an array of covariance difference eigenvalues 
                                             % computed by assuming no association of the spike train
                                             % and the stimulus.
                                             % the eigenvalues of dC (edC) will be compared to the 
                                             % distribution of null covariance eigenvalues to dermine 
                                             % which are significant 
                                             
        for ir = 1:rep                       % loop over number of repeats of shuffled spike trains 
        
            Rrnd = R(mod(randi(T)+(1:T),T)+1) ;  % random spike train obtained by shifting the real spike train by a random amounto
            iRrnd = Rrnd == 1 ;                  % with periodic boundary conditions. This preserves the burst structure contained in 
                                                 % in the VPM spike trains.
                                             
            Crnd = 1/mean(R)*cov(S1.*repmat(sqrt(Rrnd),[1 p])-repmat(mean(S1.*repmat(sqrt(Rrnd),[1 p])),[T 1])) ; 
                                                 % the null STC matrix is computed in the same way Cs is but conditioned on the 
                                                 % random spike train rather than the real one
                                             
            edCnull((ir-1)*p+(1:p)) = eig(Crnd-Cp) ;
                                                 % p null eigenvalues are computed in each iteration, so in sum
                                                 % the distributino of null eigenvalues will be constructed from 
                                                 % rep x p eigenvalues.
        end
    
        max_enull = prctile(edCnull(:),100-a) ;  % upper bound of null distribution  
        min_enull = prctile(edCnull(:),a)  ;     % lower bound of null distribution
        vprct = linspace(a,100-a,p) ;
        prctenull = prctile(edCnull(:),vprct) ;  % percentiles of null distribution (all) 
    
        isig = [find(edC>max_enull) ; find(edC<min_enull)] ;
                                                 % the significant STC filters corrspond to eigenvalues of dC outside of the null distribution 
        nsig = length(isig) ;                    % number of significant STC features
        if nsig>0 
            [~,iis] = sort(abs(edC(isig)),'descend') ;
            isig = isig(iis) ;                   % orders the STC features according to the absolute value of their corresponding eigenvalues   
        end
    
        save(['Retina_cell_' num2str(icell) '_stcsig_' stim_length{iL} '.mat'],'vdC','edC','edCnull','max_enull','min_enull','isig','nsig','a','vprct','prctenull') ;
    end
end