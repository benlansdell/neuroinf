%%%%%%%%%%%%%%%%%%%
%1 Preprocess data%
%%%%%%%%%%%%%%%%%%%

goodunits = [4,14,15,36,41];

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
nRep = 50;                      %no. sim repetitions
std = 0;
[proc, proc_withheld] = preprocess_movementinit(datafile, binsize, dt, frames, std);    
[proc, proc_withheld] = remove_bad_units(goodunits, proc, proc_withheld);
nB = size(proc.stim, 1);
fn_out = './results_coupled_GLM_L1/';
trim = 1;
pca = 0;
Dt = 20;
maxit = 20;
dt_glm = 0.1;
offset = 1;
mkdir(fn_out);

%%%%%%%%%%%%%%%%%%%%%%
%Fit L1 network model%
%%%%%%%%%%%%%%%%%%%%%%

%Run fitting...
ggs_cpl = {};
maxiter = 20;
lambdas = [.1 .3 1 3 10 30 100 300];
for l = 1:length(lambdas)
    lambda = lambdas(l);
    for icell = 1:nU
        disp(num2str(icell));
        stim = proc.stim;
        stim = stim/p;
        resp = proc.spikes{icell};
        sptrain = proc.spiketrain(:,icell);
        nicell = [(1:icell-1), (icell+1:nU)];
        %Add coupling to the other spike trains
        coupled = proc.spikes(nicell);
        for idx = 1:length(coupled)
            coupled{idx} = coupled{idx}';
        end
        stacked = proc.stacked;
        stacked = stacked/p;
        sta = stacked'*sptrain/sum(sptrain)-mean(stacked,1)'; 
        sta = reshape(sta,nF,[]);
        nspk(icell) = sum(sptrain);
        gg0 = makeFittingStruct_GLM_monkey_gauss_basisvec_refract(sta,dt,Dt);
        gg0.ihbas2 = gg0.ihbas;
        gg0.tsp = resp';
        gg0.tspi = 1;
        %Other spike trains
        gg0.tsp2 = coupled;
        %Add terms for other spike filters
        gg0.ih = zeros(size(gg0.ih,1),nU);
        opts = {'display', 'iter', 'maxiter', maxiter};
        [gg, negloglival] = MLfit_GLM_trim_L1(gg0,stim,opts,proc,trim, pca, offset, lambda);
        ggs_cpl{l,icell} = gg;
        save([fn_out '/GLM_coupled_cell_' num2str(icell) '_lambda_' num2str(lambda) '.mat'],...
            'gg'); %, 'Rt_glm');
    end
end
%Save all
save([fn_out '/all_units_network.mat'], 'ggs_cpl', 'lambdas');

for l = 1:length(lambdas)
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
    end
end

logl_glm = [];
nB = size(proc_withheld.intrial,1);
for l = 1:length(lambdas)
    for i = 1:nU
        Rt = proc_withheld.spiketrain(1:20000,i);
        tidx = ((1:nB)<=20000) & (proc_withheld.intrial==1)';
        Rt = Rt(tidx);
        %Only look within trial times...
        Rt_glm{l,i} = Rt_glm{l,i}'/nRep + 1e-8;
        if size(Rt_glm{l,i},1)==1
            Rt_glm{l,i} = Rt_glm{l,i}';
        end
        %Compute log-likelihood:
        logl_glm(l,i) = mean(Rt.*log(Rt_glm{l,i}(tidx))-(Rt_glm{l,i}(tidx))*(1/RefreshRate)) ;
    end
end

logl_glm = [];
nB = size(proc_withheld.intrial,1);
for l = 1:length(lambdas)
    for i = 1:nU
        Rt = proc_withheld.spiketrain(1:20000,i);
        %Compute log-likelihood:
        logl_glm(l,i) = mean(Rt.*log(Rt_glm{l,i})-(Rt_glm{l,i})*(1/RefreshRate)) ;
    end
end

save([fn_out '/GLM_coupled_simulation.mat'], 'Rt_glm', 'logl_glm');