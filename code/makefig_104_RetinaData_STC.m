clear ; 
cd('/Users/Shared/MANUSCRIPTS/Adrienne_Yonitan/Retinal_data') %DK files
%cd('/home/aljadeff/Documents/MATLAB/neuroinfo_retina/data') %YA files

stim_length = {'short','long'} ;
sr = 30 ;
dt = 1/sr ;
fntsz = 10 ;
for icell =  1:53

    for iL = 1:2
        load(['Retina_cell_' num2str(icell) '_stc_' stim_length{iL} '.mat']) ;
        load(['Retina_cell_' num2str(icell) '_stcsig_' stim_length{iL} '.mat']) ;
        load(['RetinaCellParameters_' stim_length{iL} '.mat']) ;
        
        p = length(edC) ;
        
        figure ; 
        
        subplot(2,1,1) ;
        ylm = 1.2*max(abs([edC(:) ; edCnull(:)]))*[-1 1] ;
        fill([0 0 p p 0],[min_enull max_enull max_enull min_enull min_enull],0.8*ones(1,3),'EdgeColor','none') ; hold on ;
        plot(setdiff(1:p,isig),edC(setdiff(1:p,isig)),'ok','MarkerFaceColor','k','MarkerSize',4) ; hold on ;
        plot(isig,edC(isig),'or','MarkerFaceColor','r','MarkerSize',4) ; hold on ;
        plot(1:p,prctenull,'Color',0.9*ones(1,3),'LineWidth',1) ;
        ylim(ylm) ;
        xlim([0 p]) ;
        set(gca, 'Layer', 'top','FontSize',fntsz) ;
        ylabel('\lambda(\Delta C)','FontSize',fntsz) ;
        xlabel('eigenvalue #','FontSize',fntsz) ;
        axis square ;
        title(num2str(icell),'FontSize',fntsz) ;
    
        subplot(2,1,2) ;
        switch kmodel
            case 1 
                plot(ctrs,Pf,'k','LineWidth',2) ; hold on ; 
                plot(ctrs,Pfs,'r','LineWidth',2) ; hold on ;
                plot(ctrs,stc_model_rect_norm(ctr)*max(Pf)/max(Psf),'b','LineWidth',2) ;
                axis tight square ; box off 
                ylim([0 1]) ;
            case 2
                imagesc(xplot,xplot,stc_model_plot,[0 1]) ; 
                xlabel('S*stc_1','FontSize',fntsz) ;
                ylabel('S*sta','FontSize',fntsz) ;

                xlim(ctr(end)*[-1 1]) ; 
                ylim(ctr(end)*[-1 1]) ; 
                axis square xy
                colormap('cool') ;
        
            case 3
                [x1,x2,x3] = ndgrid(xplot,xplot,xplot) ;
                iso02 = isosurface(x1,x2,x3,stc_model_plot,0.2) ;
                iso06 = isosurface(x1,x2,x3,stc_model_plot,0.6) ;

                p02 = patch(iso02); hold on ;
                p06 = patch(iso06); hold on ;
                % set(p02,'FaceColor',[1 0 0],'EdgeColor',[0.5 0.5 0.5],'FaceAlpha',0.3,'EdgeAlpha',0.07);
                % set(p06,'FaceColor',[0 0 1],'EdgeColor',[0.5 0.5 0.5],'FaceAlpha',0.3,'EdgeAlpha',0.07);
                set(p02,'FaceColor',[1 0 0],'EdgeColor','none','FaceAlpha',0.3);
                set(p06,'FaceColor',[0 0 1],'EdgeColor','none','FaceAlpha',0.3);

                box on ; axis square  ; 
                xlabel('S*sta','FontSize',fntsz) ;
                ylabel('S*stc_1','FontSize',fntsz) ;
                zlabel('S*stc_2','FontSize',fntsz) ;
                view([20 -40]) ;
        end
        set(gca,'FontSize',fntsz) ;

        figure ; 
        
        cm = max(abs(f(:))) ;
        for k = 1:kmodel
            subplot(kmodel,1,k) ;
            imagesc(1:Nv(icell)*pT,1:Nv(icell),reshape(f(:,k), [Nv(icell) Nv(icell)*pT]), [-cm cm]); hold on
            plot(0.5+[0 0 pT*Nv(icell) pT*Nv(icell) 0 ],0.5+[0 Nv(icell) Nv(icell) 0 0],'k') ; hold on ;
            for i = 1:pT-1
                plot(i*Nv(icell)*[1 1]+0.5,0.5+[0 Nv(icell)],'k') ; hold on ;
            end
            colormap hot;
            axis image xy off ;
            if k == 1
                title(num2str(icell),'FontSize',fntsz) ;
                ylabel('STA','FontSize',fntsz) ;
            else
                ylabel(['STC_' num2str(k-1)],'FontSize',fntsz) ;
            end
        end
     end
end  
