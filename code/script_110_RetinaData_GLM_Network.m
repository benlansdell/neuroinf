clear ; 
%cd('/Users/Shared/MANUSCRIPTS/Adrienne_Yonitan/Retinal_data') %DK files
%cd('/home/aljadeff/Documents/MATLAB/neuroinfo_retina/data') %YA files
cd('/home/lansdell/projects/neuroinf/data') %BL files
% Stimulus refresh rate (Stim frames per second)
global RefreshRate;
RefreshRate = 30; 
%Delta t of GLM, in relation to 1/RefreshRate
dt = 1;
nU = 53;
%A vector with the number of spikes for each cell 
nspk = zeros(1,nU);
logl_glm = zeros(1,nU);
rep = 10;
load(['Retina_stim_resp_network.mat']);
for icell = 1:nU
    cIR = {};
    disp(num2str(icell));% counter as loop is long
    %Prepare spikes
    nicell = [1:(icell-1), (icell+1):nU];
    resp = R(:,icell);
    coupling = R(:,nicell);
    nspk(icell) = sum(resp);
    iR = find(resp>0)';
    for idx = 1:length(nicell)
        ciR{idx} = find(coupling(:,idx)>0); 
    end
    %Prepare stim
    S = S/p;
    St = St/p;
    sta = S'*resp/sum(resp) - mean(S,1)'; 
    sta = reshape(sta,pX,pT)';
    S1 = S(:,pX*(pT-1)+1:pX*pT); 
    S1t = St(:,pX*(pT-1)+1:pX*pT); 
    %Make initial param struct   
    gg0 = makeFittingStruct_GLM_Retina(sta,min(pX,pT),dt);
    gg0.tsp = iR;
    gg0.tspi = 1;
    gg0.tsp2 = ciR;
    gg0.ih = zeros(size(gg0.ih,1),nU);
    %Compute logli of initial params
    [logli0,rr0,tt] = neglogli_GLM(gg0,S1); 
    opts = {'display', 'iter', 'maxiter', 100};
    %Do ML (requires optimization toolbox)
    [gg, negloglival] = MLfit_GLMbi(gg0,S1,opts); Rt_glm = zeros(1,Tt);
    for ir = 1:rep
        [iR_glm, vmem,Ispk] = simGLM(gg, S1t);
        Rt_glm(ceil(iR_glm)) = Rt_glm(ceil(iR_glm))+1;
    end
    Rt_glm = Rt_glm'/rep + 1e-8;
    save(['Retina_cell_' num2str(icell) '_glm_network.mat'],'gg','Rt_glm'); 
    logl_glm(icell) = mean(Rt.*log(Rt_glm)-(Rt_glm)*dt);
end
save('Retina_GLM_Network_model_All_LogLikelihood.mat','logl_glm');