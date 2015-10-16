function sim_coupled_GLM_L1(id, nRep, l)
    %Set to working directory
    wd = './';
    
    %%%%%%%%%%%%%%%%%%%
    %1 Preprocess data%
    %%%%%%%%%%%%%%%%%%%
    
    goodunits = [4,7,14,15,17,20,24,36,41];
    global RefreshRate;
    RefreshRate = 100;              %Stimulus refresh rate
    ds = 0.001;                     %Spike time resolution
    dt = ds*RefreshRate;            %Time scale spike times are resovled at,
                                    % relative to stim timescale
    datafile = './mabel_reaching_5-4-10.mat';
    nS = 4;                         %no. stim components
    nU = length(goodunits);
    frames = 80;                    %no. stim frames 
    nF = 2*frames+1;
    p = nF*nS;                      %no. stim parameters 
    binsize = 1/RefreshRate;
    standardize = 0;
    %nRep = 249;
    [proc, proc_withheld] = preprocess(datafile, binsize, dt, frames, standardize, goodunits);    
    nB = size(proc.stim, 1);
    fn_out = 'resultsL1refit/';
    trim = 1;
    Dt = 20;
    maxit = 20;
    dt_glm = 0.1;
    offset = 1;
    mkdir([wd fn_out]);
    method = 'spg';
    
    %Load fits (ggs_cpl, lambdas)
    load([fn_out '/all_units_network_L1_method_' method '.mat']);
    rng('shuffle')
    time_limit = 80;
    stim = proc_withheld.stim(:,:);
    stim = stim/p;
    
    for i = 1:nU
        ggs_cpl{l,i}.ihbas2 = ggs_cpl{l,i}.ihbas;
    end
    simstruct = makeSimStruct_GLMcpl(ggs_cpl{l,:});
    %Simulation with test stim
    disp(num2str(icell));
    Tt = size(stim,1);
    for i = 1:nU
        Rt_glm{l,i} = zeros(1,Tt);
    end
    for ir = 1:nRep
        ir
        [iR_glm,vmem,Ispk] = simGLM_monkey(simstruct, stim, time_limit);
        for i = 1:nU
            Rt_glm{l,i}(ceil(iR_glm{i})) = Rt_glm{l,i}(ceil(iR_glm{i}))+1;
        end
        save([wd fn_out '/GLM_coupled_simulation_L1_method_' method '_lambda_' num2str(lambdas(l)) '_ID_' id '.mat'], 'Rt_glm', 'nRep', 'ir');
    end
end