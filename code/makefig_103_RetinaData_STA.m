clear;
cd('/Users/Shared/MANUSCRIPTS/Adrienne_Yonitan/Retinal_data') %DK files
%cd('/home/aljadeff/Documents/MATLAB/neuroinfo_retina/data') %YA files

stim_length = {'short','long'} ;
sr = 30 ;
dt = 1/sr ;
fntsz = 10 ;
for icell = 1:53
    figure ; 
        
    for iL = 1:2
        load(['Retina_cell_' num2str(icell) '_sta_' stim_length{iL} '.mat']) ;
        load(['RetinaCellParameters_' stim_length{iL} '.mat']) ;
        
        subplot(2,2,1+2*(iL-1)) ;
    
        plot(ctrs,Pf,'k','LineWidth',2) ; hold on ; 
        plot(ctrs,Pfs,'r','LineWidth',2) ; hold on ;
        plot(ctrs,sta_model_rect_norm(ctrs)*max(Pf)/max(Psf),'b','LineWidth',2) ;
        set(gca,'FontSize',fntsz) ; axis tight square ; box off 
        ylim([0 1]) ;
        if iL==1 
            title(num2str(icell),'FontSize',fntsz) ;
        end 
        
        subplot(2,2,2+2*(iL-1)) ;
        cm = max(abs(sta)) ;
        imagesc(1:Nv(icell)*pT,1:Nv(icell),reshape(sta, [Nv(icell) Nv(icell)*pT]), [-cm cm]); hold on
    
        plot(0.5+[0 0 pT*Nv(icell) pT*Nv(icell) 0 ],0.5+[0 Nv(icell) Nv(icell) 0 0],'k') ; hold on ;
        for i = 1:pT-1
            plot(i*Nv(icell)*[1 1]+0.5,0.5+[0 Nv(icell)],'k') ; hold on ;
        end
        colormap hot;
        axis image xy off ;
     end
end    