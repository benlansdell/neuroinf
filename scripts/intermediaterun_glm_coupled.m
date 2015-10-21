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
fn_out = '/results_stampede/';
trim = 0;
Dt = 20;
maxit = 20;
dt_glm = 0.1;
mkdir([wd fn_out]);

%For each lambda value...
%Run fitting...
nU = 9;
logl_glm = [];
%Load data from the other runs and add it to Rt_glm...
r1 = load([wd fn_out 'GLM_coupled_simulation_ID_1.mat']);
r2 = load([wd fn_out 'GLM_coupled_simulation_ID_2.mat']);
r3 = load([wd fn_out 'GLM_coupled_simulation_ID_3.mat']);
r4 = load([wd fn_out 'GLM_coupled_simulation_ID_4.mat']);
r5 = load([wd fn_out 'GLM_coupled_simulation_ID_5.mat']);
r6 = load([wd fn_out 'GLM_coupled_simulation_ID_6.mat']);
r7 = load([wd fn_out 'GLM_coupled_simulation_ID_7.mat']);
r8 = load([wd fn_out 'GLM_coupled_simulation_ID_8.mat']);
r9 = load([wd fn_out 'GLM_coupled_simulation_ID_9.mat']);
r10 = load([wd fn_out 'GLM_coupled_simulation_ID_10.mat']);
r11 = load([wd fn_out 'GLM_coupled_simulation_ID_11.mat']);
r12 = load([wd fn_out 'GLM_coupled_simulation_ID_12.mat']);
r13 = load([wd fn_out 'GLM_coupled_simulation_ID_13.mat']);
r14 = load([wd fn_out 'GLM_coupled_simulation_ID_14.mat']);
r15 = load([wd fn_out 'GLM_coupled_simulation_ID_15.mat']);
r16 = load([wd fn_out 'GLM_coupled_simulation_ID_16.mat']);
r17 = load([wd fn_out 'GLM_coupled_simulation_ID_17.mat']);
r18 = load([wd fn_out 'GLM_coupled_simulation_ID_18.mat']);
r19 = load([wd fn_out 'GLM_coupled_simulation_ID_19.mat']);

%Check if Rt_glm needs normalizing
nRep = 747;
Rt_glm = {};
for idx = 1:nU
    Rt_glm{idx} = sum(r1.Rt_glm{idx},1);
    Rt_glm{idx} = Rt_glm{idx}+sum(r2.Rt_glm{idx},1);
    Rt_glm{idx} = Rt_glm{idx}+sum(r3.Rt_glm{idx},1);
    Rt_glm{idx} = Rt_glm{idx}+sum(r4.Rt_glm{idx},1);
    Rt_glm{idx} = Rt_glm{idx}+sum(r5.Rt_glm{idx},1);
    Rt_glm{idx} = Rt_glm{idx}+sum(r6.Rt_glm{idx},1);
    Rt_glm{idx} = Rt_glm{idx}+sum(r7.Rt_glm{idx},1);
    Rt_glm{idx} = Rt_glm{idx}+sum(r8.Rt_glm{idx},1);
    Rt_glm{idx} = Rt_glm{idx}+sum(r9.Rt_glm{idx},1);
    Rt_glm{idx} = Rt_glm{idx}+sum(r10.Rt_glm{idx},1);
    Rt_glm{idx} = Rt_glm{idx}+sum(r11.Rt_glm{idx},1);
    Rt_glm{idx} = Rt_glm{idx}+sum(r12.Rt_glm{idx},1);
    Rt_glm{idx} = Rt_glm{idx}+sum(r13.Rt_glm{idx},1);
    Rt_glm{idx} = Rt_glm{idx}+sum(r14.Rt_glm{idx},1);
    Rt_glm{idx} = Rt_glm{idx}+sum(r15.Rt_glm{idx},1);
    Rt_glm{idx} = Rt_glm{idx}+sum(r16.Rt_glm{idx},1);
    Rt_glm{idx} = Rt_glm{idx}+sum(r17.Rt_glm{idx},1);
    Rt_glm{idx} = Rt_glm{idx}+sum(r18.Rt_glm{idx},1);
    Rt_glm{idx} = Rt_glm{idx}+sum(r19.Rt_glm{idx},1);
end

for i = 1:nU
    Rt = proc_withheld.spiketrain(:,i);
    Rt_glm{i} = Rt_glm{i}'/nRep + 1e-8;
    if size(Rt_glm{i},1)==1
        Rt_glm{i} = Rt_glm{i}';
    end
    %Compute log-likelihood:
    logl_glm(i) = mean(Rt.*log(Rt_glm{i})-(Rt_glm{i})*(1/RefreshRate)) ;
end   
[proc, proc_withheld] = preprocess(datafile, binsize, dt, frames, standardize);    
save([wd fn_out '/preprocessed_networkglm_sims.mat'], 'proc_withheld', 'nU', 'Rt_glm', 'goodunits', 'RefreshRate')

%Save intermediate results 
save([wd fn_out '/GLM_coupled_simulation.mat'], 'Rt_glm', 'logl_glm');
load([wd fn_out '/GLM_coupled_simulation.mat'])

%Compare likelihood to uncoupled likelihood:
logl_glm_uncoupled = [];
for idx = 1:nU
    icell = goodunits(idx);
    uncoupled = load([wd '/results/GLM_cell_simulation_' num2str(icell) '.mat']);
    logl_glm_uncoupled(idx) = uncoupled.logl_glm;
end

clf
hold on 
plot(logl_glm_uncoupled, logl_glm, '.')
plot([-.5 0], [-.5 0], 'r')
xlabel('Uncoupled log-likelihood');
ylabel('Coupled log-likelihood')
saveplot(gcf, [wd fn_out '/GLM_loglikelihood_compare.eps'])