clear ;
%cd('/Users/Shared/MANUSCRIPTS/Adrienne_Yonitan/Retinal_data') %DK files
%cd('/home/aljadeff/Documents/MATLAB/neuroinfo_retina/data') %YA files
cd('/home/lansdell/projects/neuroinf/data') %BL

stim_length = {'short','long'} ;

for icell = 1:53
    for iL = 1:2
        load(['RetinaCellParameters_' stim_length{iL} '.mat']) ;
        lag = lagshifts(icell) ;
        fsize = Nv(icell)^2;
        fid = fopen(['/whitenoisec' int2str(icell) '.isk'], 'rt');
        R = textscan(fid, '%u\n');
        fclose(fid);
        R = double(R{1,1});
        T = length(R);
        fid = fopen('whitenoise.raw', 'rb');
        S = ReadFramev2(fid,T,Nx,Nv(icell),cx,x0(icell),y0(icell));
        fclose(fid);
        pX = length(S)/T;
        S = reshape(S, [pX T])';
        S = S - repmat(mean(S), [T 1]);
        S = S./repmat(std(S), [T 1]);
        R = circshift(R,-lag);
        R = R(1:end-lag);
        T = length(R);
        S = S(1:end-lag,:);
        if pT>1
            T = T - (pT-1);
            p = pX*pT;
            S1 = zeros(T,p);
            for i=1:pT
                S1(:,pX*(i-1)+1:pX*i) = S(i:T+i-1,:) ;
            end
            S = S1;
            clear S1
            R = R(pT:length(R));
        end
        itest = round(0.8*T):1:T ;
        St = S(itest,:) ;
        Rt = R(itest) ;
        S(itest,:) = [] ;
        R(itest) = [] ;
        Tt = length(itest) ;
        T = T-Tt ;
        
        save(['Retina_cell_' num2str(icell) '_stim_resp_' stim_length{iL} '.mat'],'p','pX','pT','T','Tt','S','R','St','Rt') ;
    end
end
