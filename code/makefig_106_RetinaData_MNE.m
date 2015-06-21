clear ;
cd('/Users/Shared/MANUSCRIPTS/Adrienne_Yonitan/Retinal_data') %DK files
%cd('/home/aljadeff/Documents/MATLAB/neuroinfo_retina/data') %YA files
stim_length = {'short','long'} ;

sr = 30 ;          % (Hz) sampling rate
dt = 1/sr ;       % delta t of stimulus

fntsz = 15 ;
%%
for icell = 1:53
    for iL = 1:2

        load(['Retina_cell_' num2str(icell) '_mne_' stim_length{iL} '.mat']) ;
        load(['Retina_cell_' num2str(icell) '_sta_' stim_length{iL} '.mat']) ;
        load(['Retina_cell_' num2str(icell) '_stim_resp_' stim_length{iL} '.mat']) ;
        
        ylm = 1.2*max(abs([eJ(:) ; eJnull(:)]))*[-1 1] ;
        
        figure ;
        
        subplot(nsig+1,2,[1 3 5]) ;
        fill([0 0 p p 0],[min_enull max_enull max_enull min_enull min_enull],0.8*ones(1,3),'EdgeColor','none') ; hold on ;
        plot(setdiff(1:p,isig),eJ(setdiff(1:p,isig)),'ok','MarkerFaceColor','k','MarkerSize',4) ; hold on ;
        plot(isig,eJ(isig),'or','MarkerFaceColor','r','MarkerSize',4) ; hold on ;
        ylim(ylm) ;
        xlim([0 p]) ;
        set(gca, 'Layer', 'top','FontSize',fntsz) ;
        title([num2str(icell) ', ' stim_length{iL}],'FontSize',fntsz) ;
        ylabel('\lambda(\Delta C)','FontSize',fntsz) ;
        xlabel('eigenvalue #','FontSize',fntsz) ;
        axis square ;
        
        
        subplot(nsig+1,2,2) ;
        cm = max(abs(H)) ;
        imagesc(1:sqrt(pX)*pT,1:sqrt(pX),reshape(H, [sqrt(pX) sqrt(pX)*pT]), [-cm cm]); hold on
        plot(0.5+[0 0 pT*sqrt(pX) pT*sqrt(pX) 0 ],0.5+[0 sqrt(pX) sqrt(pX) 0 0],'k') ; hold on ;
        for i = 1:pT-1
            plot(i*sqrt(pX)*[1 1]+0.5,0.5+[0 sqrt(pX)],'k') ; hold on ;
        end
        colormap hot;
        axis image xy off ;
    
   
        for j =  1:nsig
            subplot(nsig+1,2,2+2*j) ;
            
            cm = max(abs(vJ(:))) ;
            imagesc(1:sqrt(pX)*pT,1:sqrt(pX),reshape(vJ(:,isig(j)), [sqrt(pX) sqrt(pX)*pT]), [-cm cm]) ; hold on
            plot(0.5+[0 0 pT*sqrt(pX) pT*sqrt(pX) 0 ],0.5+[0 sqrt(pX) sqrt(pX) 0 0],'k') ; hold on ;
            for i = 1:pT-1
                plot(i*sqrt(pX)*[1 1]+0.5,0.5+[0 sqrt(pX)],'k') ; hold on ;
            end
            title(num2str(eJ(isig(j)),2),'FontSize',fntsz) ;
            colormap hot;
            axis image xy off ;
        end
    end
end
