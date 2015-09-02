clear ; 
stim_length = {'short','long'};
% Stimulus refresh rate (Stim frames per second)
global RefreshRate;
RefreshRate = 30; 
%Delta t of GLM, in relation to 1/RefreshRate
dt = 1/30;
%A vector with the number of spikes for each cell 
nspk = zeros(1,53);
logl_glm = zeros(2,53);
rep = 20;
trim = 1;
pcas = 0;
maxiter = 20;
offset = 1;
time_limit = 20;

%Make a fake processed structure that contains 'trial times', as a test that the trialtrim code works as it should...
processed.dt = dt;
nB = 109713;
tstarts = 1:5000:nB;
tends = tstarts + 4000;
processed.trialstartend = [tstarts', tends'];
dt_glm = 1/30;
maxit = 100;
pca = 0;
nRep = 20;

fn_out = './results3_gauss_move_5hz_maxit20_1frame_refract_testretina/';
mkdir(fn_out);

%%%%%%%%%%%%%%%%%%%%%%%%%
%2 Fitting uncoupled GLM%
%%%%%%%%%%%%%%%%%%%%%%%%%

ggs = {};
for icell = 1:10
    for iL = 1:2;
        disp(num2str(icell)) ; disp(num2str(iL)) ; % counter as loop is long
        load(['./retinadata/Retina_cell_' num2str(icell) '_stim_resp_' stim_length{iL} '.mat']) ;
        %Prepare stim
        S = S/p ;
        St = St/p ;
        sta = S'*R/sum(R) - mean(S,1)' ; 
        sta = reshape(sta,pX,pT)' ;
        S1 = S(:,pX*(pT-1)+1:pX*pT) ; 
        S1t = St(:,pX*(pT-1)+1:pX*pT) ; 
        %Prepare spikes
        nspk(icell) = sum(R) ;
        iR = find(R>0)' ;
    
        %gg0 = makeFittingStruct_GLM_monkey_gauss_basisvec_refract(sta,dt_glm,Dt);
        gg0 = makeFittingStruct_GLM_monkey(sta,dt_glm);
        gg0.tsp = iR;
        gg0.tspi = 1;
    
        %Test with random input
        %S1 = randn(size(S1));
    
        opts = {'display', 'iter', 'maxiter', maxit};
        [gg, negloglival] = MLfit_GLM_trim(gg0,S1,opts,processed,trim, pca, offset);
        ggs{icell} = gg;
    
        %Simulation with test stim
        disp(num2str(icell));
        %Test with random input
        %stim = randn(size(stim));
    
        %Simulation with test stim
        disp(num2str(icell));
        Tt = size(S1t,1);
        Rt_glm = zeros(1,Tt);
        nconverged = 0;
        for ir = 1:nRep
            ir
            [iR_glm,vmem,Ispk, conv] = simGLM_monkey(gg, S1t, time_limit);
            Rt_glm(ceil(iR_glm)) = Rt_glm(ceil(iR_glm))+1;
            nconverged = nconverged + conv;
        end
        units_conv(icell, iL) = nconverged;
        Rt_glm = Rt_glm'/nRep + 1e-8;
        save([fn_out '/Retina_cell_' num2str(icell) '_glmamended_' stim_length{iL} '.mat'], 'gg', 'Rt_glm', 'nconverged');
    end
end

%Save all
save([fn_out '/all_units.mat'], 'ggs');

%%%%%%%%%%%%%%%%%%%%%%%%%%
%3 Plot uncoupled filters%
%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_filters(ggs, proc, [fn_out '/all_units_filters.eps']);

%%%%%%%%%%%%%%%%%
%Simulate models%
%%%%%%%%%%%%%%%%%

time_limit = 40;
units_conv = zeros(nU,1);
for icell = goodunits
    load([fn_out '/GLM_cell_' num2str(icell) '.mat']);
    %Simulation with test stim
    disp(num2str(icell));
    stim = proc_withheld.stim;
    stim = stim/p;
    Tt = size(proc_withheld.stim,1);
    Rt_glm = zeros(1,Tt);
    nconverged = 0;
    for ir = 1:nRep
        ir
        [iR_glm,vmem,Ispk, conv] = simGLM_monkey(gg, stim, time_limit);
        Rt_glm(ceil(iR_glm)) = Rt_glm(ceil(iR_glm))+1;
        nconverged = nconverged + conv;
    end
    units_conv(icell) = nconverged;
    Rt_glm = Rt_glm'/nRep + 1e-8;
    save([fn_out '/GLM_cell_simulation_' num2str(icell) '.mat'], 'Rt_glm', 'nconverged');
end


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


%%%%%%%%%%%%%%%%%%%
%Fit network model%
%%%%%%%%%%%%%%%%%%%

%goodunits = [3 4 6 7 8 9 14 15 16 17 18 19 21 23 27 31 33 34 35 36 37 39 40 41];
nU = length(goodunits);
%Remove the 'bad' units from the dataset
[proc, proc_withheld] = preprocess_movementinit(datafile, binsize, dt, frames);    
[proc, proc_withheld] = remove_bad_units(goodunits, proc, proc_withheld);

%Run fitting...
ggs_cpl = {};
%For testing
maxiter = 10;
%For reals
maxiter = 100;
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
    gg0 = makeFittingStruct_GLM_monkey_gauss_basisvec(sta,dt,Dt);
    gg0.tsp = resp';
    gg0.tspi = 1;

    %Other spike trains
    gg0.tsp2 = coupled;
    %Add terms for other spike filters
    gg0.ih = zeros(size(gg0.ih,1),nU);

    opts = {'display', 'iter', 'maxiter', maxiter};
    [gg, negloglival] = MLfit_GLM_trim(gg0,stim,opts,proc,trim, pca);
    ggs_cpl{icell} = gg;

    %Simulation with test stim
    %disp(num2str(icell));
    %stim = proc_withheld.stim;
    %stim = stim/p;
    %Tt = size(proc_withheld.stim,1);
    %Rt_glm = zeros(1,Tt);
    %for ir = 1:nRep
    %    ir
    %    [iR_glm,vmem,Ispk] = simGLM_monkey(ggs{icell}, stim);
    %    Rt_glm(ceil(iR_glm)) = Rt_glm(ceil(iR_glm))+1;
    %end
    %Rt_glm = Rt_glm'/nRep + 1e-8;
    save([fn_out '/GLM_coupled_cell_' num2str(icell) '.mat'],...
        'gg'); %, 'Rt_glm');
end

%Save all
save([fn_out '/all_units_network.mat'], 'ggs_cpl');

%Plot
plot_filters_network(ggs_cpl, proc, [fn_out '/all_units_network_filters.eps']);


load([fn_out '/all_units_network.mat'])
nU = length(ggs_cpl);
%Simulate network model...
stim = proc_withheld.stim;
stim = stim/p;
clear proc;
clear proc_withheld;
for icell = 1:nU
    %Simulation with test stim
    disp(num2str(icell));
    Tt = size(stim,1);
    Rt_glm = zeros(1,Tt);
    for ir = 1:nRep
        ir
        [iR_glm,vmem,Ispk] = simGLM_monkey(ggs_cpl{icell}, stim);
        Rt_glm(ceil(iR_glm)) = Rt_glm(ceil(iR_glm))+1;
    end
    Rt_glm = Rt_glm'/nRep + 1e-8;
    save([fn_out '/GLM_coupled_simulation_' num2str(icell) '.mat'], 'Rt_glm');
end