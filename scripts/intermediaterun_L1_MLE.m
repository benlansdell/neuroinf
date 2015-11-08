%Set to working directory
wd = '.';

%%%%%%%%%%%%%%%%%%%
%1 Preprocess data%
%%%%%%%%%%%%%%%%%%%
goodunits = [4,7,14,15,17,20,24,36,41];
lambdas = [.1 .3 1 3 10 30 100 300];

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
fn_out = '/results_L1_MLE/';
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
for l = 3:length(lambdas)
    l
    lambda = lambdas(l);
    r1 = load([wd fn_out 'GLM_coupled_simulation_L1_method_spg_lambda_' num2str(lambda) '_ID_1.mat']);
    r2 = load([wd fn_out 'GLM_coupled_simulation_L1_method_spg_lambda_' num2str(lambda) '_ID_2.mat']);
    r3 = load([wd fn_out 'GLM_coupled_simulation_L1_method_spg_lambda_' num2str(lambda) '_ID_3.mat']);
    %Check if Rt_glm needs normalizing
    Rt_glm = {};
    for idx = 1:nU
        Rt_glm{idx} = sum(r1.Rt_glm{l,idx},1);
        Rt_glm{idx} = Rt_glm{idx}+sum(r2.Rt_glm{l,idx},1);
        Rt_glm{idx} = Rt_glm{idx}+sum(r3.Rt_glm{l,idx},1);
    end
    nB = size(Rt_glm{1}, 2);

    for i = 1:nU
        Rt = proc_withheld.spiketrain(1:nB,i);
        Rt_glm{i} = Rt_glm{i}'/nRep + 1e-8;
        if size(Rt_glm{i},1)==1
            Rt_glm{i} = Rt_glm{i}';
        end
        %Compute log-likelihood:
        logl_glm(l, i) = mean(Rt.*log(Rt_glm{i})-(Rt_glm{i})*(1/RefreshRate)) ;
    end   
    %[proc, proc_withheld] = preprocess(datafile, binsize, dt, frames, standardize);   
    %proc_withheld 
    save([wd fn_out '/preprocessed_networkglm_sims_lambda_' num2str(lambda) '.mat'], 'proc_withheld', 'nU', 'Rt_glm', 'goodunits', 'RefreshRate')
    coh_out = ['coherence_lambda_' num2str(lambda)];
    fn_ins{l-2} = ['/preprocessed_networkglm_sims_lambda_' num2str(lambda) '.mat'];
    jackknifecoherence([wd fn_out], ['/preprocessed_networkglm_sims_lambda_' num2str(lambda) '.mat'], coh_out)
end

coh_out_all = ['fullsim'];
jackknifecoherence_all([wd fn_out], fn_ins, coh_out_all)

%Save intermediate results 
%save([wd fn_out '/GLM_coupled_simulation.mat'], 'Rt_glm', 'logl_glm');
%load([wd fn_out '/GLM_coupled_simulation.mat'])

%Compare likelihood to uncoupled likelihood:
logl_glm_uncoupled = [];
for idx = 1:nU
    icell = goodunits(idx);
    uncoupled = load([wd '/results/GLM_cell_simulation_' num2str(icell) '.mat']);
    logl_glm_uncoupled(idx) = uncoupled.logl_glm;
end

clf
semilogx(lambdas(3:end), sum(logl_glm(3:end,:),2))
xlabel('\lambda');
ylabel('Total coupled log-likelihood')
title('Full simulation')
saveplot(gcf, [wd fn_out '/GLM_loglikelihood_compare_semilog.eps'])

maxlogl = max(logl_glm(3:end-1,:), [], 1);
uncoupledlogl = logl_glm_uncoupled;
clf
plot(uncoupledlogl, maxlogl, 'o')
hold on
plot([-.4 0], [-.4 0], 'r')
title('Full simulation')
xlabel('Uncoupled log-likelihood');
ylabel('Coupled log-likelihood')
saveplot(gcf, [wd fn_out '/GLM_loglikelihood_compare.eps'])

logl = logl_glm(3:end,:);
clf
cmap = jet(7);
hold all
for idx = 1:6
    plot(1:size(logl,2), logl(idx,:), 'o', 'Color', cmap(idx,:))
end
plot(1:size(logl,2), uncoupledlogl, 'o', 'Color', cmap(7,:))
legend('\lambda = 1 (fully coupled)', '\lambda = 3', '\lambda = 10', '\lambda = 30', '\lambda = 100', '\lambda = 300', 'uncoupled', 'location', 'eastoutside')
xlabel('Unit');
ylabel('Coupled log-likelihood')
title('Full simulation')
saveplot(gcf, [wd fn_out '/GLM_loglikelihood_lambda.eps'], 'eps', [7 4])


logl = logl_glm(3:end,:);
clf
cmap = jet(6);
hold all
for idx = 1:6
    plot(1:size(logl,2), logl(idx,:)-uncoupledlogl, 'o', 'Color', cmap(idx,:))
end
plot([0,10],[0 0], 'k--')
legend('\lambda = 1 (fully coupled)', '\lambda = 3', '\lambda = 10', '\lambda = 30', '\lambda = 100', '\lambda = 300', 'location', 'eastoutside')
xlabel('Unit');
ylabel('Coupled - uncoupled log-likelihood')
title('Full simulation')
saveplot(gcf, [wd fn_out '/GLM_loglikelihood_lambda_relative.eps'], 'eps', [7 4])
