clear ; 
%cd('/Users/Shared/MANUSCRIPTS/Adrienne_Yonitan/Retinal_data') %DK files
%cd('/home/aljadeff/Documents/MATLAB/neuroinfo_retina/data') %YA files
cd('/home/lansdell/projects/neuroinf/data') %BL files

stim_length = {'short','long'} ;

dt = 1/30 ;%ds/sr ;    % delta t of stimulus 
sr = 1/dt ; 
dt = 1;
nspk = zeros(1,53) ;       % a vector with the number of spikes for each cell 

logl_glm = zeros(2,53) ;

% Filter_rank = 5 ; % Number of column/row vector pairs to use

rep = 10 ;
    
global RefreshRate;  % Stimulus refresh rate (Stim frames per second)
RefreshRate = sr ; 

%for icell = 1:53
%    for iL = 1:2
icell = 1; iL = 1;
        disp(num2str(icell)) ; disp(num2str(iL)) ; % counter as loop is long

        load(['Retina_cell_' num2str(icell) '_stim_resp_' stim_length{iL} '.mat']) ;
        S = S/p ;
        St = St/p ;
        sta = S'*R/sum(R) - mean(S,1)' ; 
        sta = reshape(sta,pX,pT)' ;
        %S1 = S(:,1:pX) ; 
        %S1t = St(:,1:pX) ;
        S1 = S(:,pX*(pT-1)+1:pX*pT) ; 
        S1t = St(:,pX*(pT-1)+1:pX*pT) ; 
        nspk(icell) = sum(R) ;
        iR = find(R>0)' ;
    
        % Compute STA and use as initial guess for k
    
        % 4. Do ML fitting of params with simulated data %=====================

        %  Initialize params for fitting --------------
        gg0 = makeFittingStruct_GLM_Retina(sta,min(pX,pT),dt);
        gg0.tsp = iR ;
        gg0.tspi = 1 ;
        [logli0,rr0,tt] = neglogli_GLM(gg0,S1); % Compute logli of initial params

        % Do ML estimation of model params
        opts = {'display', 'iter', 'maxiter', 100};
        [gg, negloglival] = MLfit_GLMbi(gg0,S1,opts); % do ML (requires optimization toolbox)
 
        Rt_glm = zeros(1,Tt) ;
    
        for ir = 1:rep
            [iR_glm, vmem,Ispk] = simGLM(gg, S1t) ;
            Rt_glm(ceil(iR_glm)) = Rt_glm(ceil(iR_glm))+1 ;
        end
        Rt_glm = Rt_glm'/rep + 1e-8 ;

        save(['Retina_cell_' num2str(icell) '_glmtest_' stim_length{iL} '.mat'],'gg','Rt_glm') ; 

        logl_glm(iL,icell) = mean(Rt.*log(Rt_glm)-(Rt_glm)*dt) ;
%    end
%end
 save('Retina_GLMtest_model_All_LogLikelihood.mat','logl_glm') ;