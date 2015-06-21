clear ; 
cd('/Users/Shared/MANUSCRIPTS/Adrienne_Yonitan/Retinal_data') %DK files
%cd('/home/aljadeff/Documents/MATLAB/neuroinfo_retina/data') %YA files

stim_length = {'short','long'} ;

dt = 1/30 ;        % delta t of stimulus 

nspk = zeros(1,53) ;       % a vector with the number of spikes for each cell 

rep = 500 ;        % number of repeats used to compute the null distribution of covariance eigenvalues

a = 0.2 ;          % level of significance required to call an STC feature significant
                   % if an eigenvalue of dC is smaller than the a-th percentile of the null eigenvalue 
                   % distribution or larger than the 100-a percentile, it will be determined to be significant

nbns = 10 ;        % number of bins for each relevant stimulus dimension
                   % the stimulus frames (and correspoding spikes) will be 
                   % binned into nbins^nsig bins 
nplot = nbns*3 ;  

linearrect = @(x) x.*(sign(x)+1)/2 ;

for icell = 1:53; 
    for iL = 1:2
        disp(num2str(icell)) ; disp(num2str(iL)) ; % counter as loop is long

        load(['Retina_cell_' num2str(icell) '_stim_resp_' stim_length{iL} '.mat']) ;
        load(['Retina_cell_' num2str(icell) '_sta_' stim_length{iL} '.mat']) ;
        load(['Retina_cell_' num2str(icell) '_stcsig_' stim_length{iL} '.mat']) ;
        
        iR = R>0 ;
        kmodel = 1 + min(nsig,2) ; 
        
        f = [sta vdC(:,isig(1:kmodel-1))] ;
        
        z = S*f ;
        zt = St*f ;
        zx = max(abs(z(:)))+1e-3 ;
        bns = linspace(-zx,zx,nbns) ;
        xplot = linspace(-zx,zx,nplot) ;
        dbns = bns(2)-bns(1) ;
        ctr = bns(1:nbns-1)+dbns/2 ;

        switch kmodel
            case 1 
                Pf  = histcn(z,bns) ;
                Pfs = histcn(z(iR,:),bns) ;
            case 2 
                Pf  = histcn(z,bns,bns) ;
                Pfs = histcn(z(iR,:),bns,bns) ;
            case 3
                Pf  = histcn(z,bns,bns,bns) ;
                Pfs = histcn(z(iR,:),bns,bns,bns) ;
        end
        % the STC feature is the eigenvector of dC with largest eigenvalue (in absolute value)
        
        Ps = mean(iR) ;                      % Ps  - probability of a spike
                                             %     - the number of spikes divided by the number of stimulus 'frames'
        
        Pf = Pf/(T*dbns^kmodel) ;            %     - a normalized count of the projections of stimuli on the STA, STC1 in nbns^2 bins
        Pfs = Pfs/(T*Ps*dbns^kmodel) ;    %     - a normalized count of the projections of stimuli on STA, STC1 in nbns^2 bins conditioned on a spike
        
        Psf = Pfs./Pf*Ps ;                   % Psf - probability of a spike given stimulus projected on a feature
                                             %     - computed by invoking Bayes' rule (see details in text)
        Psf(isnan(Psf)) = 0 ;                % setting NaNs to 0
        
        switch kmodel
            case 1
                stc_model = @(x)interpn(ctr,Psf,x,'linear') ;
                x1 = xplot ;
                X = x1 ;
            case 2 
                stc_model = @(x)interpn(ctr,ctr,Psf,x(:,1),x(:,2),'linear') ;
                [x1,x2] = ndgrid(xplot,xplot) ;
                x1v = reshape(x1,nplot^kmodel,1) ;
                x2v = reshape(x2,nplot^kmodel,1) ;
                X = [x1v x2v] ;
            case 3
                stc_model = @(x)interpn(ctr,ctr,ctr,Psf,x(:,1),x(:,2),x(:,3),'linear') ;
                [x1,x2,x3] = ndgrid(xplot,xplot,xplot) ;
                x1v = reshape(x1,nplot^kmodel,1) ;
                x2v = reshape(x2,nplot^kmodel,1) ;
                x3v = reshape(x2,nplot^kmodel,1) ;
                X = [x1v x2v x3v] ;
        end
                                             % first the STC model is defined to be a 
                                             % smoothing spline interpolation of the spiking
                                             % nonlinearity computed using binned probability distributions  
        stc_model_rect = @(x) linearrect(stc_model(x))+1e-8 ;
                                             % then the model passed through a linear rectifier, to 
                                             % assure that the probability of a spike is non-negative
                                             % the infinitesimal constant offset (1e-8) ensures that the 
                                             % loglikelihood calculation below does not diverge
        stc_model_z = stc_model_rect(z) ;
        stc_model_z(isnan(stc_model_z)) = 0 ;
        Ps_model = sum(stc_model_z) ;        % finally the model has to be renormalized such that the number of 
                                             % of predicted spikes (on the training set) is fixed to what was 
                                             % found in the experiment. This step is needed because of the 
                                             % smoothing and rectifying operations  
        stc_model_rect_norm = @(x) stc_model_rect(x)*Ps*T/Ps_model ;
                                             % this normalization is done by computing the number of predicted 
                                             % spikes from the smooth rectified model (Ps_model)                                                                                                 
                                             % then the model multiplied by the ratio of the actual number of spikes and Ps_model
        stc_model_plot = reshape(stc_model_rect_norm(X),nplot*ones(1,kmodel)) ;
        
        save(['Retina_cell_' num2str(icell) '_stc_' stim_length{iL} '.mat'],'stc_model_rect_norm','f','kmodel','Pfs','Psf','Pf','stc_model_plot','bns','ctr','xplot') ; 
        
        %logl_stc(i) = mean(Rt.*log(stc_model_rect_norm(zt))-stc_model_rect_norm(zt)*(ds/sr)) ;
    end
end