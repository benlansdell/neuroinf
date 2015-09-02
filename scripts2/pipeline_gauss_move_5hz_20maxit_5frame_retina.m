clear ; 
stim_length = {'short','long'};
% Stimulus refresh rate (Stim frames per second)
global RefreshRate;
RefreshRate = 30; 
%Delta t of GLM, in relation to 1/RefreshRate
dt = 1/30;
%A vector with the number of spikes for each cell 
nspk = zeros(1,53);
logl_glm = zeros(2,53);
rep = 10;
trim = 1;
pcas = 0;
maxiter = 20;


%Make a fake processed structure that contains 'trial times', as a test that the trialtrim code works as it should...
processed.dt = dt;
nB = 109713;
tstarts = 1:5000:nB;
tends = tstarts + 4000;
processed.trialstartend = [tstarts', tends'];

for icell = 1:53
    for iL = 1:2
        disp(num2str(icell)) ; disp(num2str(iL)) ; % counter as loop is long
        load(['./retinadata/Retina_cell_' num2str(icell) '_stim_resp_' stim_length{iL} '.mat']) ;
        %Prepare stim
        S = S/p ;
        St = St/p ;
        sta = S'*R/sum(R) - mean(S,1)' ; 
        sta = reshape(sta,pX,pT)' ;
        S1 = S(:,pX*(pT-1)+1:pX*pT) ; 
        S1t = St(:,pX*(pT-1)+1:pX*pT) ; 
        %Prepare spikes
        nspk(icell) = sum(R) ;
        iR = find(R>0)' ;
        %Make initial param struct   
        gg0 = makeFittingStruct_GLM_retina_verify(sta,dt);
        gg0.tsp = iR ;
        gg0.tspi = 1 ;
        %Compute logli of initial params
        [logli0,rr0,tt] = neglogli_GLM(gg0,S1); 
        opts = {'display', 'iter', 'maxiter', maxiter};
        %Do ML (requires optimization toolbox)
        [gg, negloglival] = MLfit_GLM_trim(gg0,S1,opts, processed, trim, pcas);
        Rt_glm = zeros(1,Tt) ;
        for ir = 1:rep
            [iR_glm, vmem,Ispk] = simGLM_monkey(gg, S1t) ;
            Rt_glm(ceil(iR_glm)) = Rt_glm(ceil(iR_glm))+1 ;
        end
        Rt_glm = Rt_glm'/rep + 1e-8 ;
        save(['./retinadata/Retina_cell_' num2str(icell) '_glmamended_' stim_length{iL} '.mat'],'gg','Rt_glm') ; 
        logl_glm(iL,icell) = mean(Rt.*log(Rt_glm)-(Rt_glm)*dt) ;
    end
end
save('./retinadata/Retina_GLMamedned_model_All_LogLikelihood.mat','logl_glm') ;