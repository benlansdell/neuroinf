clear ;
%cd('/Users/Shared/MANUSCRIPTS/Adrienne_Yonitan/Retinal_data') %DK files
%cd('/home/aljadeff/Documents/MATLAB/neuroinfo_retina/data') %YA files
cd('/home/lansdell/projects/neuroinf/data') %BL files

pltpath = '/home/lansdell/projects/neuroinf/plots/';
plotwidth = 8;
plotheight = 8;

stim_length = {'short','long'} ;

sr = 30 ;          % (Hz) sampling rate
dt = 1/sr ;       % delta t of spike train

fntsz = 15 ;

for icell = 1:53
    for iL = 1:2
        fn_out = [pltpath '/Retina_cell_' num2str(icell) '_glm_' stim_length{iL} '.eps'];
        load(['Retina_cell_' num2str(icell) '_glm_' stim_length{iL} '.mat']) ;
        load(['Retina_cell_' num2str(icell) '_sta_' stim_length{iL} '.mat']) ;
        
        figure ;

        subplot(2,1,1) ;
        cm = max(abs(gg.k(:))) ;
        Splot = zeros(sqrt(pX),sqrt(pX)*pT) ;
        for i = 1:pT
            Splot(1:sqrt(pX),(1:sqrt(pX))+(i-1)*sqrt(pX)) = reshape(gg.k(i,:),sqrt(pX),sqrt(pX)) ;
        end
        imagesc(1:sqrt(pX)*pT,1:sqrt(pX),Splot, [-cm cm]) ; hold on
        plot(0.5+[0 0 pT*sqrt(pX) pT*sqrt(pX) 0 ],0.5+[0 sqrt(pX) sqrt(pX) 0 0],'k') ; hold on ;
        for i = 1:pT-1
            plot(i*sqrt(pX)*[1 1]+0.5,0.5+[0 sqrt(pX)],'k') ; hold on ;
        end
        colormap hot;
        axis image xy off ;
        set(gca,'FontSize',fntsz,'Yaxislocation','right') ;
        
        subplot(2,1,2) ;
        plot(gg.iht,gg.ihbas*gg.ih)
        plot(dt*gg.iht,zeros(size(gg.iht)),'k') ; hold on ;
        plot(dt*gg.iht,gg.ihbas*gg.ih,'b','LineWidth',3) ; axis tight ; box off ;
        legend('GLM spike history filter') ;
        xlabel('time relative to spike (s)') ;
        set(gca,'FontSize',fntsz) ;
        saveplot(gcf, fn_out, 'eps', [plotwidth plotheight]);
    end
end