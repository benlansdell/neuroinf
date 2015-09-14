%%%%%%%%%%%%%%%%%%%
%1 Preprocess data%
%%%%%%%%%%%%%%%%%%%
goodunits = [4,7,8,14,15,16,17,20,24,36,37,41];

%no. frames:
nframes = [40 60 80 100 120];

global RefreshRate;
RefreshRate = 100;              %Stimulus refresh rate
ds = 0.001;                     %Spike time resolution
dt = ds*RefreshRate;            %Time scale spike times are resovled at,
                                % relative to stim timescale
datafile = './mabel_reaching_5-4-10.mat';
nU = 45;                        %no. units
nS = 4;                         %no. stim components
frames = 40;                   %no. stim frames 
nF = 2*frames+1;
p = nF*nS;                      %no. stim parameters 
binsize = 1/RefreshRate;
nRep = 20;                      %no. sim repetitions
std = 0;
fn_out = './results5_gauss_move_5hz_maxit20_BIC/';
trim = 1;
pca = 0;
Dt = 20;
maxit = 20;
dt_glm = 0.1;
mkdir(fn_out);

%%%%%%%%%%%%%%%%%%%%%%%%%
%2 Fitting uncoupled GLM%
%%%%%%%%%%%%%%%%%%%%%%%%%
for frames = nframes
    nF = 2*frames+1;
    [proc, proc_withheld] = preprocess_movementinit(datafile, binsize, dt, frames, std);    
    nB = size(proc.stim, 1);
    
    ggs = {};
    for icell = goodunits
        disp(num2str(icell));
        stim = proc.stim;
        stim = stim/p;
        resp = proc.spikes{icell};
        sptrain = proc.spiketrain(:,icell);
    
        stacked = proc.stacked;
        stacked = stacked/p;
        sta = stacked'*sptrain/sum(sptrain)-mean(stacked,1)'; 
        sta = reshape(sta,nF,[]);
    
        nspk(icell) = sum(sptrain);
        gg0 = makeFittingStruct_GLM_monkey_gauss_basisvec_refract(sta,dt_glm,Dt);
        gg0.tsp = resp';
        gg0.tspi = 1;
    
        opts = {'display', 'iter', 'maxiter', maxit};
        [gg, negloglival] = MLfit_GLM_trim(gg0,stim,opts,proc,trim, pca);
        ggs{icell} = gg;
        plot_filters(ggs, proc, [fn_out '/unit_' num2str(icell) '_filters.eps'], [icell]);
        save([fn_out '/GLM_cell_' num2str(icell) '_frames_' num2str(frames) '.mat'],...
            'gg', 'negloglival'); %, 'Rt_glm');
    end
    
    %Save all
    save([fn_out '/all_units_frames_' num2str(frames) '.mat'], 'ggs');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %3 Plot uncoupled filters%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_filters(ggs, proc, [fn_out '/all_units_filters_frames_' num2str(frames) '.eps'], goodunits);
    
    %%%%%%%%%%%%%%%%%
    %Simulate models%
    %%%%%%%%%%%%%%%%%
    
    %time_limit = 40;
    %units_conv = zeros(nU,1);
    %for icell = goodunits
    %    load([fn_out '/GLM_cell_' num2str(icell) '.mat']);
    %    %Simulation with test stim
    %    disp(num2str(icell));
    %    stim = proc_withheld.stim;
    %    stim = stim/p;
    %    %stim = stim*0;
    %    Tt = size(proc_withheld.stim,1);
    %    Rt_glm = zeros(1,Tt);
    %    nconverged = 0;
    %    for ir = 1:nRep
    %        ir
    %        [iR_glm,vmem,Ispk, conv] = simGLM_monkey(gg, stim, time_limit);
    %        Rt_glm(ceil(iR_glm)) = Rt_glm(ceil(iR_glm))+1;
    %        nconverged = nconverged + conv;
    %    end
    %    units_conv(icell) = nconverged;
    %    Rt_glm = Rt_glm'/nRep + 1e-8;
    %    save([fn_out '/GLM_cell_simulation_' num2str(icell) '.mat'], 'Rt_glm', 'nconverged');
    %end
end

%Plot likelihood curves
lhood = zeros(length(nframes), length(goodunits));
frames = 40;
[proc, proc_withheld] = preprocess_movementinit(datafile, binsize, dt, frames, std);
k = (1+5+[1 3 5 7 9])*nU;
n = size(proc.stim, 1);
for i = 1:length(nframes)
    frames = nframes(i);
    for j = 1:length(goodunits)
        icell = goodunits(j);
        load([fn_out '/GLM_cell_' num2str(icell) '_frames_' num2str(frames) '.mat']);
        lhood(i,j) = negloglival;
    end
end
sumlhood = sum(lhood,2);
BIC = 2*sumlhood + k'*log(n);
AIC = 2*sumlhood + 2*k';

clf
%plot(nframes, sum(lhood,2))
ncomp = [1 3 5 7 9];
plot(ncomp, BIC, ncomp, AIC)
xlabel('no. stim components')
ylabel('negative log likelihood')
title('Acausal+causal stim. filter (stim <> spikes)')
legend('BIC', 'AIC')
saveplot(gcf, [fn_out 'lhood.eps'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Comparison to other filter types%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lhood_acausal = zeros(length(nframes), length(goodunits));
for i = 1:length(nframes)
    frames = nframes(i);
    for j = 1:length(goodunits)
        icell = goodunits(j);
        load(['./results5_gauss_move_5hz_maxit20_BIC_nooffset_acausal/GLM_cell_' num2str(icell) '_frames_' num2str(frames) '.mat']);
        lhood_acausal(i,j) = negloglival;
    end
end

lhood_causal = zeros(length(nframes), length(goodunits));
for i = 1:length(nframes)
    frames = nframes(i);
    for j = 1:length(goodunits)
        icell = goodunits(j);
        load(['./results5_gauss_move_5hz_maxit20_BIC_nooffset/GLM_cell_' num2str(icell) '_frames_' num2str(frames) '.mat']);
        lhood_causal(i,j) = negloglival;
    end
end

clf
nU = length(goodunits);
plot(1:nU, -lhood(3,:), 'r', 1:nU, -lhood_causal(3,:), 'b', 1:nU, -lhood_acausal(3,:), 'k');
legend('Both', 'Causal', 'Acausal')
xlabel('Unit')
ax = gca;
set(ax, 'XTick', 1:12);
set(ax,'XTickLabel',num2str(goodunits'));
ylabel('log-likelihood')
saveplot(gcf, [fn_out '/lhood_comparison.eps'])

clf
nU = length(goodunits);
plot(1:nU, -lhood(3,:)+lhood_causal(3,:), 'b', 1:nU, -lhood_acausal(3,:)+lhood_causal(3,:), 'k');
legend('Both-Causal', 'Acausal-Causal')
xlabel('Unit')
ax = gca;
set(ax, 'XTick', 1:12);
set(ax,'XTickLabel',num2str(goodunits'));
ylabel('\Delta log-likelihood')
saveplot(gcf, [fn_out '/lhood_comparison_relative.eps'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4 Plot uncoupled simulations%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for icell = goodunits
    clf
    sigma_fr = .25;
    sigma_fr = sigma_fr*RefreshRate;
    sz = sigma_fr*3*2;
    x = linspace(-sz/2, sz/2, sz);
    gaussFilter_fr = exp(-x.^2/(2*sigma_fr^2));
    gaussFilter_fr = gaussFilter_fr/sum(gaussFilter_fr);
    
    tstart = ceil(600*RefreshRate);
    tend = ceil(750*RefreshRate);
    fn_in = [fn_out '/GLM_cell_' num2str(icell) '.mat'];
    load([fn_out '/GLM_cell_' num2str(icell) '.mat'])
    load([fn_out '/GLM_cell_simulation_' num2str(icell) '.mat'])
    tidx = tstart:tend;
    truesp = proc_withheld.spiketrain(tidx,icell);
    pillowsimsp = Rt_glm(tidx);
    subplot(2,1,1)
    hold on
    plot(tidx*proc.binsize, 500*proc.grip(tidx), 'k');
    plot(tidx*proc.binsize, proc.cursor(tidx,1), 'b');
    plot(tidx*proc.binsize, proc.cursor(tidx,2), 'Color', [0 0.6 0]);
    plot(tidx*proc.binsize, proc.cursor(tidx,3), 'r');
    legend('Grip', 'Curs x', 'Curs Y', 'Curs Z')
    subplot(2,1,2)
    gftruesp = conv(truesp, gaussFilter_fr, 'same');
    gfpillowsimsp = conv(pillowsimsp, gaussFilter_fr, 'same');
    plot(tidx*proc.binsize, gftruesp, tidx*proc.binsize, gfpillowsimsp);
    title('n rep: 20')
    xlabel('seconds')
    ylabel('predicted probability spiking')
    legend('true spike train', 'Pillow''s GLM')
    saveplot(gcf, [fn_out '/GLM_cell_' num2str(icell) '_sim.eps'], 'eps', [6 6]);
    
    clf
    sigma_fr = .01;
    sigma_fr = sigma_fr*RefreshRate;
    sz = sigma_fr*3*2;
    x = linspace(-sz/2, sz/2, sz);
    gaussFilter_fr = exp(-x.^2/(2*sigma_fr^2));
    gaussFilter_fr = gaussFilter_fr/sum(gaussFilter_fr);
    
    tstart = ceil(600*RefreshRate);
    tend = ceil(605*RefreshRate);
    tidx = tstart:tend;
    truesp = proc_withheld.spiketrain(tidx,icell);
    pillowsimsp = Rt_glm(tidx);
    subplot(2,1,1)
    hold on
    plot(tidx*proc.binsize, 500*proc.grip(tidx), 'k');
    plot(tidx*proc.binsize, proc.cursor(tidx,1), 'b');
    plot(tidx*proc.binsize, proc.cursor(tidx,2), 'Color', [0 0.6 0]);
    plot(tidx*proc.binsize, proc.cursor(tidx,3), 'r');
    legend('Grip', 'Curs x', 'Curs Y', 'Curs Z')
    subplot(2,1,2)
    gftruesp = conv(truesp, gaussFilter_fr, 'same');
    gfpillowsimsp = conv(pillowsimsp, gaussFilter_fr, 'same');
    spidx = truesp==1;
    plot(tidx(spidx)*proc.binsize, truesp(spidx)-.95, '.', tidx*proc.binsize, gfpillowsimsp);
    title('n rep: 20')
    xlabel('seconds')
    ylabel('predicted probability spiking')
    legend('true spike train', 'Pillow''s GLM')
    saveplot(gcf, [fn_out '/GLM_cell_' num2str(icell) '_sim_zoom.eps'], 'eps', [6 6]);
end
