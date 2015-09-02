function makefig_RetinaData_GLM(pltpath)    
    plotwidth = 8;
    plotheight = 8;
    stim_length = {'short','long'} ;
    sr = 30 ;          % (Hz) sampling rate
    dt = 1/sr ;       % delta t of spike train
    fntsz = 15 ;
    for icell = 1:6
        for iL = 1:2
            fn_out = [pltpath '/Retina_cell_' num2str(icell) '_glmamended_' stim_length{iL} '.eps'];
            load(['./retinadata/Retina_cell_' num2str(icell) '_sta_' stim_length{iL} '.mat']) ;
            load(['./retinadata/Retina_cell_' num2str(icell) '_stim_resp_' stim_length{iL} '.mat']) ;
            yonatanglm = load(['~/projects/neuroinf_retina_monkey_2015_08_14/data/Retina_cell_' num2str(icell) '_glm_' stim_length{iL} '.mat']) ;
            yonatanglm_amend = load(['~/projects/neuroinf_retina_monkey_2015_08_14/data/Retina_cell_' num2str(icell) '_glmamended_' stim_length{iL} '.mat']) ;

            try
                ownglm = load([pltpath '/Retina_cell_' num2str(icell) '_glmamended_' stim_length{iL} '.mat']) ;
            catch
                display(['Can''t find ' pltpath '/Retina_cell_' num2str(icell) '_glmamended_' stim_length{iL} '.mat. Continuing'])
                continue 
            end
            figure;    
            subplot(2,1,1);
            cm = max(abs(ownglm.gg.k(:)));
            Splot = zeros(sqrt(pX),sqrt(pX)*pT);
            for i = 1:pT
                Splot(1:sqrt(pX),(1:sqrt(pX))+(i-1)*sqrt(pX)) = reshape(ownglm.gg.k(i,:),sqrt(pX),sqrt(pX));
            end
            imagesc(1:sqrt(pX)*pT,1:sqrt(pX),Splot, [-cm cm]); hold on
            plot(0.5+[0 0 pT*sqrt(pX) pT*sqrt(pX) 0 ],0.5+[0 sqrt(pX) sqrt(pX) 0 0],'k') ; hold on ;
            for i = 1:pT-1
                plot(i*sqrt(pX)*[1 1]+0.5,0.5+[0 sqrt(pX)],'k') ; hold on ;
            end
            colormap hot;
            axis image xy off ;
            set(gca,'FontSize',fntsz,'Yaxislocation','right') ;
            
            subplot(2,1,2) ;
            plot(ownglm.gg.iht,ownglm.gg.ihbas*ownglm.gg.ih)
            plot(dt*ownglm.gg.iht,zeros(size(ownglm.gg.iht)),'k') ; hold on ;
            plot(dt*ownglm.gg.iht,ownglm.gg.ihbas*ownglm.gg.ih,'b','LineWidth',3) ; axis tight ; box off ;
            legend('GLM spike history filter') ;
            xlabel('time relative to spike (s)') ;
            set(gca,'FontSize',fntsz) ;
            saveplot(gcf, fn_out, 'eps', [plotwidth plotheight]);
            figure 
    
            %Plot simulation results
            RefreshRate = 30;
            sigma_fr = .1;
            sigma_fr = sigma_fr*RefreshRate;
            sz = sigma_fr*3*2;
            x = linspace(-sz/2, sz/2, sz);
            gaussFilter_fr = exp(-x.^2/(2*sigma_fr^2));
            gaussFilter_fr = gaussFilter_fr/sum(gaussFilter_fr);
        
            tstart = 1501;
            tend = 1750;
            ii = tstart:tend;
            tt = ii*dt;

            gftruesp = conv(Rt, gaussFilter_fr, 'same');
            gfownsp = conv(ownglm.Rt_glm, gaussFilter_fr, 'same');
            gfyonatonamendsimsp = conv(yonatanglm_amend.Rt_glm, gaussFilter_fr, 'same');
            gfyonatonsimsp = conv(yonatanglm.Rt_glm, gaussFilter_fr, 'same');

            plot(tt, Rt(ii), tt, gfownsp(ii),tt, gfyonatonamendsimsp(ii),tt, gfyonatonsimsp(ii));
            legend('true spikes', 'my glm', 'yonatan glm', 'yonatan glm amended')
            xlabel('seconds')
            ylabel('prob of spiking')
            saveplot(gcf, [fn_out(1:end-4) '_sim.eps']);        
        end
    end
end