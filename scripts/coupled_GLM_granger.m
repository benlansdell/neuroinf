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
fn_out = './results_coupled_GLM_granger/';
trim = 1;
pca = 0;
Dt = 20;
maxit = 20;
dt_glm = 0.1;
offset = 1;
mkdir(fn_out);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fit Granger network model%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Run fitting on entire model
ggs_cpl = {};
maxiter = 20;
grangerstat = zeros(nU, nU);

for icell = 1:nU
    disp(num2str(icell));
    resp = proc.spikes{icell};
    sptrain = proc.spiketrain(:,icell);

    stim = proc.stim;
    stim = stim/p;
    stacked = proc.stacked;
    stacked = stacked/p;
    sta = stacked'*sptrain/sum(sptrain)-mean(stacked,1)'; 
    sta = reshape(sta,nF,[]);
    nspk(icell) = sum(sptrain);
    gg0 = makeFittingStruct_GLM_monkey_gauss_basisvec_refract(sta,dt,Dt);
    gg0.ihbas2 = gg0.ihbas;
    gg0.tsp = resp';
    gg0.tspi = 1;

    %Add coupling to the other spike trains
    nicell = [(1:icell-1), (icell+1:nU)];
    coupled = proc.spikes(nicell);
    for idx = 1:length(coupled)
        coupled{idx} = coupled{idx}';
    end
    gg0.tsp2 = coupled;
    gg0.ih = zeros(size(gg0.ih,1),nU);

    %Run optimization
    opts = {'display', 'iter', 'maxiter', maxiter};
    [gg, negloglival_all] = MLfit_GLM_trim(gg0,stim,opts,proc,trim, pca, offset);
    ggs_cpl{icell} = gg;
    save([fn_out '/GLM_coupled_cell_' num2str(icell) '.mat'], 'gg');

    %Run by systematically dropping coupling terms from fitting
    nC = nU-1;
    for j = 1:nC
        jcell = nicell(j);
        nijcell = [(1:j-1), (j+1:nC)];
        coupled = proc.spikes(nicell(nijcell))';
        for idx = 1:length(coupled)
            coupled{idx} = coupled{idx}';
        end
        gg0.tsp2 = coupled;
        gg0.ih = zeros(size(gg0.ih,1),nC);
        opts = {'display', 'iter', 'maxiter', maxiter};
        [gg, negloglival] = MLfit_GLM_trim(gg0,stim,opts,proc,trim, pca, offset);
        grangerstat(icell, jcell) = negloglival - negloglival_all;
    end
end

%Perform hypothesis test, 