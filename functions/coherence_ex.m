clear ;
cd('/Users/Shared/MANUSCRIPTS/Adrienne_Yonitan/VPM_data') ; % location of .mat files
% modified on 16 June by DK 
% modified on 15 June by YA

c_sta   = [  0   0 128]/255 ;
c_pwsta = [128 128 255]/255 ;
c_stc   = [  0 128   0]/255 ;
c_mne   = [128   0 255]/255 ;
c_glm   = [255 128   0]/255 ;
c_tun   = [128 128 128]/255 ;
white   = [255 255 255]/255 ;
c_stm   = [0 0 0] ;
LL = zeros(7,6) ;

Names = {'37_A2' ; '46_BC' ; '57_C4' ; '83_E2' ; '88_E1' ; '92_D2' ; '93_C4' } ;
lgnd = {'37 A2' ; '46 BC' ; '57 C4' ; '83 E2' ; '88 E1' ; '92 D2' ; '93 C4' } ;

sr = 1000 ;        % (Hz) sampling rate
p = 150 ;          % stimulus dimensionality 
ds = 2 ;           % downsampling factor: one in every ds (=2) whisker position data points will be included in the analysis 
dt = ds/sr ;       % delta t of stimulus 

fmax=25;            % maximum plot frequency (set for convenience)


for i = 1:7 ; %set to cell 57C4 - use 1:7 for all cells
    
    % Predicted spike train from XXX model
    load(['VPM_cell_' Names{i} '_stim_resp.mat']) ;
    St = St(:,1) ;

    load(['VPM_cell_' Names{i} '_PredictedSpikeTrains.mat']) ;
    
    nT = length(Rt) ;
    t = dt*(1:Tt) ; % just the time base
    
    K=0.0; % just a switch to look at LL compared to null hypothesis
    
    LL(i,7) = sum(Rt.*log(mean(Rt)*dt)-mean(Rt)*dt-log(factorial(Rt)))  ;
    LL(i,1) = sum(Rt.*log(Rt_sta.*dt)-Rt_sta.*dt-log(factorial(Rt))) - K*LL(i,7) ;
    LL(i,2) = sum(Rt.*log(Rt_stak.*dt)-Rt_stak.*dt-log(factorial(Rt))) - K*LL(i,7) ;
    LL(i,3) = sum(Rt.*log(Rt_stc.*dt)-Rt_stc.*dt-log(factorial(Rt))) - K*LL(i,7)  ;
    LL(i,4) = sum(Rt.*log(Rt_mne.*dt)-Rt_mne.*dt-log(factorial(Rt))) - K*LL(i,7) ;
    LL(i,5) = sum(Rt.*log(Rt_glm.*dt)-Rt_glm.*dt-log(factorial(Rt))) - K*LL(i,7)  ;
    LL(i,6) = sum(Rt.*log(Rt_tun'.*dt)-Rt_tun'.*dt-log(factorial(Rt))) - K*LL(i,7)  ;


NW=12 ; % NW = time-bandwidth product so 12/20s = 0.6 Hz ; this is a soft number!
params.tapers = [NW NW-1];
params.Fs = 1/dt ; % sampling frequency (inverse of sampling time)
params.fpass=[0 fmax]; % plot from 0 to 2 Hz
params.pad = 4; % 4-times padding, a rule of thumb
params.err=[0 0.05]; % 0.05 confidence interval
params.trialave = 0; % one trace
    

    
    [SRt, ~] = mtspectrumc(Rt-mean(Rt),params); % other Chronux routines are set up for point processes rather than continuous data
    [SRt_sta, ~] = mtspectrumc(Rt_sta-mean(Rt_sta),params);
    [SRt_stak, ~] = mtspectrumc(Rt_stak-mean(Rt_stak),params);
    [SRt_stc, ~] = mtspectrumc(Rt_stc-mean(Rt_stc),params);
    [SRt_mne, ~] = mtspectrumc(Rt_mne-mean(Rt_mne),params);
    [SRt_glm, ~] = mtspectrumc(Rt_glm-mean(Rt_glm),params);
    [SRt_tun, ~] = mtspectrumc(Rt_tun-mean(Rt_tun),params);
    [SSt, ~] = mtspectrumc(St(:,1)-mean(St(:,1)),params);
    
    SRt = SRt/norm(SRt) ;
    SRt_sta = SRt_sta/norm(SRt_sta) ;
    SRt_stak = SRt_stak/norm(SRt_stak) ;
    SRt_stc = SRt_stc/norm(SRt_stc) ;
    SRt_mne = SRt_mne/norm(SRt_mne) ;
    SRt_glm = SRt_glm/norm(SRt_glm) ;
    SRt_tun = SRt_tun/norm(SRt_tun) ;
    SSt     = SSt/norm(SSt) ;
    
NW=25 ; % NW = time-bandwidth product so 20/20s = 1.25 Hz ; this is a soft number!
params.tapers = [NW NW-1];
params.Fs = 1/dt ; % sampling frequency (inverse of sampling time)
params.fpass=[0 fmax]; % plot from 0 to 2 Hz
params.pad = 4; % 4-times padding, a rule of thumb
params.err=[2 0.02]; % 0.05 confidence interval
params.trialave = 0; % one trace

    [Coh_sta,phi_sta,~,~,~,freq,confC,~,~] = coherencyc(Rt-mean(Rt),Rt_sta-mean(Rt_sta),params); % other Chronux routines are set up for point processes rather than continuous data, or mixed point and continuous; S12, S1, S2 are cross-power and respective powers for the same NW - but We estimated these for a smaller NW; f is the frequency base
    [Coh_stak,phi_stak,~,~,~,~,~,~,~] = coherencyc(Rt-mean(Rt),Rt_stak-mean(Rt_stak),params); % other Chronux routines are set up for point processes rather than continuous data, or mixed point and continuous; S12, S1, S2 are cross-power and respective powers for the same NW - but We estimated these for a smaller NW; f is the frequency base
    [Coh_stc,phi_stc,~,~,~,~,~,~,~] = coherencyc(Rt-mean(Rt),Rt_stc-mean(Rt_stc),params); % other Chronux routines are set up for point processes rather than continuous data, or mixed point and continuous; S12, S1, S2 are cross-power and respective powers for the same NW - but We estimated these for a smaller NW; f is the frequency base
    [Coh_mne,phi_mne,~,~,~,~,~,~,~] = coherencyc(Rt-mean(Rt),Rt_mne-mean(Rt_mne),params); % other Chronux routines are set up for point processes rather than continuous data, or mixed point and continuous; S12, S1, S2 are cross-power and respective powers for the same NW - but We estimated these for a smaller NW; f is the frequency base
    [Coh_glm,phi_glm,~,~,~,~,~,~,~] = coherencyc(Rt-mean(Rt),Rt_glm-mean(Rt_glm),params); % other Chronux routines are set up for point processes rather than continuous data, or mixed point and continuous; S12, S1, S2 are cross-power and respective powers for the same NW - but We estimated these for a smaller NW; f is the frequency base
    [Coh_tun,phi_tun,~,~,~,~,~,~,~] = coherencyc(Rt-mean(Rt),Rt_tun-mean(Rt_tun),params); % other Chronux routines are set up for point processes rather than continuous data, or mixed point and continuous; S12, S1, S2 are cross-power and respective powers for the same NW - but We estimated these for a smaller NW; f is the frequency base
    
    figure ;
    
    subplot(2,1,1) ;
    plot(dt*(1:Tt),St(:,1),'k','Linewidth',1) ;
    axis tight

    subplot(2,1,2) ;
    plot(t,Rt     ,'Color',[200 200 200]/255,'LineWidth',1) ; hold on ;
    plot(t,Rt_sta ,'Color',c_sta,'LineWidth',1) ; hold on ;
    plot(t,Rt_stak,'Color',c_pwsta,'LineWidth',1) ; hold on ;
    plot(t,Rt_stc ,'Color',c_stc,'LineWidth',1) ; hold on ;
    plot(t,Rt_mne ,'Color',c_mne,'LineWidth',1) ; hold on ;
    plot(t,Rt_glm ,'Color',c_glm,'LineWidth',1) ; hold on ;
    plot(t,Rt_tun ,'Color',c_tun,'LineWidth',1) ; hold on ;
    
    xlabel('t [s]') ;
    ylabel('p(spike)') ;
    legend('data','sta','pwsta','stc','mne','glm','tuning') ;
    axis tight
    
    figure ;
    subplot(4,1,1) ;
    
    plot(freq,log(SRt),'Color',[200 200 200]/255,'LineWidth',1) ; hold on ;
    plot(freq,log(SRt_sta),'Color',c_sta,'LineWidth',1) ; hold on ;
    plot(freq,log(SRt_stak),'Color',c_pwsta,'LineWidth',1) ; hold on ;
    plot(freq,log(SRt_stc),'Color',c_stc,'LineWidth',1) ; hold on ;
    plot(freq,log(SRt_mne),'Color',c_mne,'LineWidth',1) ; hold on ;
    plot(freq,log(SRt_glm),'Color',c_glm,'LineWidth',1) ; hold on ;
    plot(freq,log(SRt_tun),'Color',c_tun,'LineWidth',1) ; hold on ;
    plot(freq,log(SSt),'Color',c_stm,'LineWidth',2) ; hold on ;
    ylabel('Log(Spectral Power)','FontSize',10) ;
   
    subplot(4,1,2);
    plot(freq,(Coh_sta),'Color',c_sta,'LineWidth',1) ; hold on ;
    plot(freq,(Coh_stak),'Color',c_pwsta,'LineWidth',1) ; hold on ;
    plot(freq,(Coh_stc),'Color',c_stc,'LineWidth',1) ; hold on ;
    plot(freq,(Coh_mne),'Color',c_mne,'LineWidth',1) ; hold on ;
    plot(freq,(Coh_glm),'Color',c_glm,'LineWidth',1) ; hold on ;
    plot(freq,(Coh_tun),'Color',c_tun,'LineWidth',1) ; hold on ;
    plot([0 fmax],confC*ones(1,2),'Color',c_sta,'LineWidth',2) ; hold on ; % plot confidence interval
    ylabel('Magnitude of Coherence - |C|','FontSize',10) ;
    ylim([0 1]) ;

    phi_sta_s = mod(phi_sta,2*pi) ; 
    phi_sta_s(phi_sta_s>pi) = phi_sta_s(phi_sta_s>pi) - 2*pi ; 
    phi_sta_s(Coh_sta<confC) = NaN ;
    
    phi_stak_s = mod(phi_stak,2*pi) ; 
    phi_stak_s(phi_stak_s>pi) = phi_stak_s(phi_stak_s>pi) - 2*pi ; 
    phi_stak_s(Coh_stak<confC) = NaN ;
    
    phi_stc_s = mod(phi_stc,2*pi) ; 
    phi_stc_s(phi_stc_s>pi) = phi_stc_s(phi_stc_s>pi) - 2*pi ; 
    phi_stc_s(Coh_stc<confC) = NaN ;
    
    phi_mne_s = mod(phi_mne,2*pi) ; 
    phi_mne_s(phi_mne_s>pi) = phi_mne_s(phi_mne_s>pi) - 2*pi ; 
    phi_mne_s(Coh_mne<confC) = NaN ;
    
    phi_glm_s = mod(phi_glm,2*pi) ; 
    phi_glm_s(phi_glm_s>pi) = phi_glm_s(phi_glm_s>pi) - 2*pi ; 
    phi_glm_s(Coh_glm<confC) = NaN ;
    
    phi_tun_s = mod(phi_tun,2*pi) ; 
    phi_tun_s(phi_tun_s>pi) = phi_tun_s(phi_tun_s>pi) - 2*pi ; 
    phi_tun_s(Coh_tun<confC) = NaN ;
    
    subplot(4,1,3); 
    
    plot(freq,phi_sta_s,'Color',c_sta,'LineWidth',1) ; hold on ;
    plot(freq,phi_stak_s,'Color',c_pwsta,'LineWidth',1) ; hold on ;
    plot(freq,phi_stc_s,'Color',c_stc,'LineWidth',1) ; hold on ;
    plot(freq,phi_mne_s,'Color',c_mne,'LineWidth',1) ; hold on ;
    plot(freq,phi_glm_s,'Color',c_glm,'LineWidth',1) ; hold on ;
    plot(freq,phi_tun_s,'Color',c_tun,'LineWidth',1) ; hold on ;
    
    ylabel('Phase of Coherence - arg{C}','FontSize',10) ;
    set(gca,'ytick',-pi:pi/2:pi,'yticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'}) ;
    ylim([-pi pi]) ;
    
    subplot(4,1,4) ;
    plot(freq,(Coh_sta).^2,'Color',c_sta,'LineWidth',1) ; hold on ;
    plot(freq,(Coh_stak).^2,'Color',c_pwsta,'LineWidth',1) ; hold on ;
    plot(freq,(Coh_stc).^2,'Color',c_stc,'LineWidth',1) ; hold on ;
    plot(freq,(Coh_mne).^2,'Color',c_mne,'LineWidth',1) ; hold on ;
    plot(freq,(Coh_glm).^2,'Color',c_glm,'LineWidth',1) ; hold on ;
    plot(freq,(Coh_tun).^2,'Color',c_tun,'LineWidth',1) ; hold on ;
    ylabel('Probability - |C|^2','FontSize',10) ;

    xlabel('frequency(Hz)','FontSize',10) ;
  
    
end
%%
figure ; 
barh = bar(1:7,LL) ; 
barh(1).FaceColor = c_sta   ;
barh(2).FaceColor = c_pwsta ;
barh(3).FaceColor = c_stc   ; 
barh(4).FaceColor = c_mne   ;
barh(5).FaceColor = c_glm   ;
barh(6).FaceColor = c_tun   ;
barh(7).FaceColor = white   ;

legend('sta','pwsta','stc','mne','glm','tuning','constant') ;
    
set(gca,'xticklabel',lgnd,'FontSize',15) ;
ylabel('loglikelihood','FontSize',15) ;
