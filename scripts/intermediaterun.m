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
nRep = 747;                      %no. sim repetitions
standardize = 0;
[proc, proc_withheld] = preprocess(datafile, binsize, dt, frames, standardize, goodunits);    
nB = size(proc.stim, 1);
fn_out = '/results_intermed747/';
trim = 0;
Dt = 20;
maxit = 20;
dt_glm = 0.1;
mkdir([wd fn_out]);

%Load data from the other runs and add it to Rt_glm...
r1 = load([wd fn_out '/GLM_coupled_simulation1.mat']);
r2 = load([wd fn_out '/GLM_coupled_simulation2.mat']);
r3 = load([wd fn_out '/GLM_coupled_simulation3.mat']);
r4 = load([wd fn_out '/GLM_coupled_simulation4.mat']);
r5 = load([wd fn_out '/GLM_coupled_simulation5.mat']);
r6 = load([wd fn_out '/GLM_coupled_simulation6.mat']);
r7 = load([wd fn_out '/GLM_coupled_simulation7.mat']);
r8 = load([wd fn_out '/GLM_coupled_simulation8.mat']);
r9 = load([wd fn_out '/GLM_coupled_simulation9.mat']);

nU = 9;
%Check if Rt_glm needs normalizing
Rt_glm = {};
for idx = 1:nU
	Rt_glm{idx} = r1.Rt_glm{idx};
	Rt_glm{idx} = Rt_glm{idx}+r2.Rt_glm{idx};
	Rt_glm{idx} = Rt_glm{idx}+r3.Rt_glm{idx};
	Rt_glm{idx} = Rt_glm{idx}+r4.Rt_glm{idx};
	Rt_glm{idx} = Rt_glm{idx}+r5.Rt_glm{idx};
    Rt_glm{idx} = Rt_glm{idx}+r6.Rt_glm{idx};
    Rt_glm{idx} = Rt_glm{idx}+r7.Rt_glm{idx};
    Rt_glm{idx} = Rt_glm{idx}+r8.Rt_glm{idx};
    Rt_glm{idx} = Rt_glm{idx}+r9.Rt_glm{idx};
end

nRep = 747;

logl_glm = [];
for i = 1:nU
    Rt = proc_withheld.spiketrain(:,i);
    Rt_glm{i} = Rt_glm{i}'/nRep + 1e-8;
    if size(Rt_glm{i},1)==1
        Rt_glm{i} = Rt_glm{i}';
    end
    %Compute log-likelihood:
    logl_glm(i) = mean(Rt.*log(Rt_glm{i})-(Rt_glm{i})*(1/RefreshRate)) ;
end

%Save intermediate results 
save([wd fn_out '/GLM_coupled_simulation.mat'], 'Rt_glm', 'logl_glm');

%Compare likelihood to uncoupled likelihood:
load([wd fn_out '/GLM_coupled_simulation.mat'])
logl_glm_uncoupled = [];
for idx = 1:nU
    icell = goodunits(idx);
    uncoupled = load([wd fn_out '/GLM_cell_simulation_' num2str(icell) '.mat']);
    logl_glm_uncoupled(idx) = uncoupled.logl_glm;
end
clf
scatter(logl_glm_uncoupled, logl_glm)
hold on; plot([-1.5 0], [-1.5 0], 'r')
xlabel('Uncoupled log-likelihood');
ylabel('Coupled log-likelihood')
saveplot(gcf, [wd fn_out '/GLM_loglikelihood_compare.eps'])

%Save for computation of coherence 
[proc, proc_withheld] = preprocess(datafile, binsize, dt, frames, standardize);    
save([wd fn_out '/preprocessed_networkglm_sims.mat'], 'proc_withheld', 'nU', 'Rt_glm', 'goodunits', 'RefreshRate')

%Compuate coherence 