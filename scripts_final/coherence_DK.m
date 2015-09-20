clear ;
%cd('/Users/Shared/MANUSCRIPTS/Adrienne_Yonitan/simulationtraces') ;% set to location of .mat files on YOUR computer

 
sr = 100 ;        % (Hz) sampling rate
fmax=5;            % maximum plot frequency (set for convenience)

 
%load(['GLM_sims_cell_4.mat']) %108237 samples
load(['./results4_gauss_move_5hz_maxit20_5frame/GLM_sims_cell_8_inmovement2.mat']) % 169 spikes
    %Chronux parameter set-up
NW=5 ; % NW = time-bandwidth product so 5/1.7 = 3 Hz ; this is a soft number!
params.tapers = [NW NW-1];
params.Fs = sr ; % sampling frequency (inverse of sampling time)
params.fpass=[0 fmax]; % plot from 0 to 2 Hz
params.pad = 4; % 4-times padding, a rule of thumb
params.err=[1 0.05]; % 0.05 confidence interval
params.trialave = 0; % one trace

 
    [S,f,Ser] = mtspectrumc(truesp-mean(truesp),params); % other Chronux routines are set up for point processes rather than continuous data
    [Sm,f,Ser] = mtspectrumc(simsp-mean(simsp),params);
    [Sn,f,Ser] = mtspectrumc(simsp_cpl-mean(simsp),params);

       
    figure ; % plot data and models here
    plot(f,log(S),'r','LineWidth',1) ; hold on ;
     plot(f,log(Sm),'g','LineWidth',1) ; hold on ;
    plot(f,log(Sn),'b','LineWidth',1) ; hold on ;
    ylabel('Log(Pwr)','FontSize',8) ;

   

        
NW=5 ; % NW = time-bandwidth product (p) so 50/1080s = 0.05 Hz ; this is a soft number!
params.tapers = [NW NW-1];
params.Fs = sr ; % sampling frequency (inverse of sampling time)
params.fpass=[0 fmax]; % plot from 0 to 2 Hz
params.pad = 4; % 4-times padding, a rule of thumb
params.err=[2 0.02]; % 0.02 confidence interval
params.trialave = 0; % one trace

 
    [Coh_m,phi_m,~,~,~,freq,confC,~,~] = coherencyc(truesp-mean(truesp),simsp-mean(simsp),params); % other Chronux routines are set up for point processes rather than continuous data, or mixed point and continuous; S12, S1, S2 are cross-power and respective powers for the same NW - but We estimated these for a smaller NW; f is the frequency base
    [Coh_n,phi_n,~,~,~,~,~,~,~] = coherencyc(truesp-mean(truesp),simsp_cpl-mean(simsp_cpl),params); % other Chronux routines are set up for point processes rather than continuous data, or mixed point and continuous; S12, S1, S2 are cross-power and respective powers for the same NW - but We estimated these for a smaller NW; f is the frequency base

    
   figure ; % plot data and models here
    plot(freq,(Coh_m),'g','LineWidth',1) ; hold on ;
    plot(freq,(Coh_n),'b','LineWidth',1) ; hold on ;
    plot([0 fmax],confC*ones(1,2),'k','LineWidth',2) ; hold on ; % plot confidence interval
    xlabel('Frequency [Hz]','FontSize',10) ;