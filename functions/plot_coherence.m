function plot_coherence(fn_out, Rt, Rt_glmcpl, Rt_glm, dt)
    c_glm   = [255 128   0]/255 ;
    c_glmcpl   = [0 128   255]/255 ;

    fmax = 250;
    NW=12 ; % NW = time-bandwidth product so 12/20s = 0.6 Hz ; this is a soft number!
    params.tapers = [NW NW-1];
    params.Fs = 1/dt ; % sampling frequency (inverse of sampling time)
    params.fpass=[0 fmax]; % plot from 0 to 2 Hz
    params.pad = 4; % 4-times padding, a rule of thumb
    params.err=[0 0.05]; % 0.05 confidence interval
    params.trialave = 0; % one trace
        
    [SRt, ~] = mtspectrumc(Rt-mean(Rt),params); % other Chronux routines are set up for point processes rather than continuous data
    [SRt_glm, ~] = mtspectrumc(Rt_glm-mean(Rt_glm),params);
    [SRt_glmcpl, ~] = mtspectrumc(Rt_glmcpl-mean(Rt_glmcpl),params);
        
    SRt = SRt/norm(SRt);
    SRt_glm = SRt_glm/norm(SRt_glm);
    SRt_glmcpl = SRt_glmcpl/norm(SRt_glmcpl);
        
    NW=25 ; % NW = time-bandwidth product so 20/20s = 1.25 Hz ; this is a soft number!
    params.tapers = [NW NW-1];
    params.Fs = 1/dt ; % sampling frequency (inverse of sampling time)
    params.fpass=[0 fmax]; % plot from 0 to 2 Hz
    params.pad = 4; % 4-times padding, a rule of thumb
    params.err=[2 0.02]; % 0.05 confidence interval
    params.trialave = 0; % one trace
    
    [Coh_glm,phi_glm,~,~,~,freq,confC,~,~] = coherencyc(Rt-mean(Rt),Rt_glm-mean(Rt_glm),params); % other Chronux routines are set up for point processes rather than continuous data, or mixed point and continuous; S12, S1, S2 are cross-power and respective powers for the same NW - but We estimated these for a smaller NW; f is the frequency base
    [Coh_glmcpl,phi_glmcpl,~,~,~,~,~,~,~] = coherencyc(Rt-mean(Rt),Rt_glmcpl-mean(Rt_glmcpl),params); % other Chronux routines are set up for point processes rather than continuous data, or mixed point and continuous; S12, S1, S2 are cross-power and respective powers for the same NW - but We estimated these for a smaller NW; f is the frequency base

    phi_glm_s = mod(phi_glm,2*pi); 
    phi_glm_s(phi_glm_s>pi) = phi_glm_s(phi_glm_s>pi) - 2*pi; 
    phi_glm_s(Coh_glm<confC) = NaN;
       
    phi_glmcpl_s = mod(phi_glmcpl,2*pi); 
    phi_glmcpl_s(phi_glmcpl_s>pi) = phi_glmcpl_s(phi_glmcpl_s>pi) - 2*pi; 
    phi_glmcpl_s(Coh_glmcpl<confC) = NaN;

    subplot(3,1,1);
    plot(freq,(Coh_glm),'Color',c_glm,'LineWidth',1); hold on;
    plot(freq,(Coh_glmcpl),'Color',c_glmcpl,'LineWidth',1);
    plot([0 fmax],confC*ones(1,2),'Color',c_glm,'LineWidth',2); % plot confidence interval
    ylabel('Magnitude of Coherence - |C|','FontSize',10);
    ylim([0 1]);
            
    subplot(3,1,2); 
    
    plot(freq,phi_glm_s,'Color',c_glm,'LineWidth',1); hold on;
    plot(freq,phi_glmcpl_s,'Color',c_glmcpl,'LineWidth',1);
    ylabel('Phase of Coherence - arg{C}','FontSize',10);
    set(gca,'ytick',-pi:pi/2:pi,'yticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'});
    ylim([-pi pi]);
    
    subplot(3,1,3);
    plot(freq,(Coh_glm).^2,'Color',c_glm,'LineWidth',1); hold on;
    plot(freq,(Coh_glmcpl).^2,'Color',c_glmcpl,'LineWidth',1);
    saveplot(gcf, [fn_out], 'eps', [5 8])
end