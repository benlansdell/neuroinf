%Set to working directory
wd = '.';

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
nU = length(goodunits);
nS = 4;                         %no. stim components
frames = 80;                    %no. stim frames 
nF = 2*frames+1;
p = nF*nS;                      %no. stim parameters 
binsize = 1/RefreshRate;
nRep = 747;                     %no. sim repetitions
standardize = 0;
[proc, proc_withheld] = preprocess(datafile, binsize, dt, frames, standardize, goodunits);    
nB = size(proc.stim, 1);
trim = 1;
Dt = 20;
maxit = 20;
dt_glm = 0.1;

ggs_cpl = {};
maxiter = 20;
icell = 1;

disp(['Fitting unit ' num2str(goodunits(icell))]);
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
gg0 = makeFittingStruct_GLM_monkey(sta,dt,Dt);
gg0.ihbas2 = gg0.ihbas;
gg0.tsp = resp';
gg0.tspi = 1;
%Other spike trains
gg0.tsp2 = coupled;
%Add terms for other spike filters
gg0.ih = zeros(size(gg0.ih,1),nU);
opts = {'display', 'iter', 'maxiter', maxiter};
gg = gg0; Stim = stim; optimArgs = opts; 
MAXSIZE  = 1e7;  % Maximum amount to be held in memory at once;
offset = 1;
opts = optimset('Gradobj','on','Hessian','on', optimArgs{:});
% Set initial params
prs0 = extractFitPrs_GLM_trim(gg,Stim,MAXSIZE,proc, trim, offset);

%Compute log likelihood (involves using spikeconv_mex)
[fv,gradval,H] = Loss_GLM_logli(prs0);
display(['Log likelihood of initial params given spikes is ' num2str(fv)])

%Returns
%
%americano: Log likelihood of initial params given spikes is 582.022 
%stampede: Log likelihood of initial params given spikes is 582.022

%Same.... so spikeconv_mex works on stampede despite being compiled on americano?
%Or it's compiled just in time?