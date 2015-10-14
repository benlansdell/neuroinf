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
fn_out = '/results_L1_stampede/';
trim = 0;
Dt = 20;
maxit = 20;
dt_glm = 0.1;
mkdir([wd fn_out]);

%For each lambda value...
%Run fitting...
nU = 9;
lambdas = [.1 .3 1 3 10 30 100 300];
logl_glm = [];
for l = 1:length(lambdas)
    lambda = lambdas(l);
    %Load data from the other runs and add it to Rt_glm...
    r1 = load([wd fn_out 'GLM_coupled_simulation_L1_method_spg_lambda_' num2str(lambda) '_ID_1.mat']);
    r2 = load([wd fn_out 'GLM_coupled_simulation_L1_method_spg_lambda_' num2str(lambda) '_ID_2.mat']);
    r3 = load([wd fn_out 'GLM_coupled_simulation_L1_method_spg_lambda_' num2str(lambda) '_ID_3.mat']);

    %Check if Rt_glm needs normalizing
    nRep = 747;
    Rt_glm = {};
    for idx = 1:nU
        Rt_glm{idx} = r1.Rt_glm{l,idx};
        Rt_glm{idx} = Rt_glm{idx}+r2.Rt_glm{l,idx};
        Rt_glm{idx} = Rt_glm{idx}+r3.Rt_glm{l,idx};
    end

    for i = 1:nU
        Rt = proc_withheld.spiketrain(:,i);
        Rt_glm{i} = Rt_glm{i}'/nRep + 1e-8;
        if size(Rt_glm{i},1)==1
            Rt_glm{i} = Rt_glm{i}';
        end
        %Compute log-likelihood:
        logl_glm(l, i) = mean(Rt.*log(Rt_glm{i})-(Rt_glm{i})*(1/RefreshRate)) ;
    end   
    [proc, proc_withheld] = preprocess(datafile, binsize, dt, frames, standardize);    
    save([wd fn_out '/preprocessed_networkglm_sims_lambda_' num2str(lambda) '.mat'], 'proc_withheld', 'nU', 'Rt_glm', 'goodunits', 'RefreshRate')
end

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
for l = 1:length(lambdas)
    plot(logl_glm(end,:), logl_glm(l,:), '.')
    plot([-2.7 0], [-2.7 0], 'r')
end
xlabel('Uncoupled log-likelihood');
ylabel('Coupled log-likelihood')
legend('.1', '.3', '1', '3', '10', '30', '100', '300');
saveplot(gcf, [wd fn_out '/GLM_loglikelihood_compare.eps'])

clf
semilogx(lambdas, sum(logl_glm,2))
xlabel('\lambda')
ylabel('Total log-likelihood')
saveplot(gcf, [wd fn_out '/GLM_loglikelihood_compare_semilog.eps'])

%Save for computation of coherence 
