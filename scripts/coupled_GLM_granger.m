%Set to working directory
wd = '.';

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
nRep = 249;                      %no. sim repetitions
standardize = 0;
[proc, proc_withheld] = preprocess(datafile, binsize, dt, frames, standardize, goodunits);    
nB = size(proc.stim, 1);
fn_out = './results/';
trim = 1;
Dt = 20;
maxit = 20;
dt_glm = 0.1;
offset = 1;
alphalevel = 0.001;
mkdir([wd fn_out]);

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
    gg0 = makeFittingStruct_GLM_monkey(sta,dt,Dt);
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
    [gg, negloglival_all] = MLfit_GLM_trim(gg0,stim,opts,proc,trim, offset);
    ggs_cpl{icell} = gg;

    %Run by exhaustively dropping coupling terms from fitting
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
        [gg, negloglival] = MLfit_GLM_trim(gg0,stim,opts,proc,trim, offset);
        grangerstat(icell, jcell) = 2*(negloglival - negloglival_all);
    end
end

%Perform hypothesis test, 
dof = size(gg0.ih,1);
GCpval = 1-chi2cdf(grangerstat, dof);
GCsig = multiple_sig(GCpval, alphalevel);

save([wd fn_out '/grangerstats.mat'], 'grangerstat', 'GCpval', 'GCsig');
