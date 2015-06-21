clear ;
cd('/Users/Shared/MANUSCRIPTS/Adrienne_Yonitan/Retinal_data') %DK files
%cd('/home/aljadeff/Documents/MATLAB/neuroinfo_retina/data') %YA files
fntsz = 10 ;
stim_length = {'short','long'} ;
for icell = 1:53
    for iL = 1:2
        
        figure ; 
    
        load(['Retina_cell_' num2str(icell) '_stim_resp_' stim_length{iL} '.mat']) ;
    
        C = cov(S) ;
        [vC,eC] = eig(C) ;
        [eC,iC] = sort(diag(eC),'descend') ;
        vC = vC(:,iC) ;

        plot(1:p,eC,'ok','MarkerSize',8,'MarkerFaceColor','w') ;
        xlabel('eigenvalue #','FontSize',fntsz) ;
        ylabel('\lambda','FontSize',fntsz) ;
        set(gca,'FontSize',fntsz) ;
        ylim([0 1.2]) ;
        
        figure ;
        for j = 1:12
            subplot(12,1,j) ;
            imagesc(reshape(vC(:,j), sqrt(pX), pT*sqrt(pX)), [-0.25 0.25]); hold on
            plot(0.5+[0 0 pT*sqrt(pX) pT*sqrt(pX) 0 ],0.5+[0 sqrt(pX) sqrt(pX) 0 0],'k') ; hold on ;
            for k = 1:pT-1
                plot(k*sqrt(pX)*[1 1],0.5+[0 sqrt(pX)],'k') ; hold on ;
            end
            axis equal off ;
        end
        colormap('gray')    
    end
end