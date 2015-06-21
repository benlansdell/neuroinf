clear ;
cd('/Users/Shared/MANUSCRIPTS/Adrienne_Yonitan/Retinal_data') %DK files
%cd('/home/aljadeff/Documents/MATLAB/neuroinfo_retina/data') %YA files

load('RetinaCellParameters_long.mat') ;
pT = 1 ; 

x0net = min(x0) ;
y0net = min(y0) ;
Nv = 10 + max(max(x0)-x0net,max(y0)-y0net) ;
lag = 0 ;
fsize = Nv^2;

fid = fopen(['/whitenoisec' int2str(1) '.isk'], 'rt');
Ri = textscan(fid, '%u\n');
fclose(fid);
Ri = double(Ri{1,1});
T = length(Ri);

R = zeros(T,53) ;
for icell = 1:53
    fid = fopen(['/whitenoisec' int2str(icell) '.isk'], 'rt');
    Ri = textscan(fid, '%u\n');
    fclose(fid);
    R(:,icell) = double(Ri{1,1});
end

fid = fopen('whitenoise.raw', 'rb');
S = ReadFramev2(fid,T,Nx,Nv,cx,x0net,y0net);
fclose(fid);

pX = length(S)/T;
p = pT*pX ;

S = reshape(S, [pX T])';
S = S - repmat(mean(S), [T 1]);
S = S./repmat(std(S), [T 1]);


itest = round(0.8*T):1:T ;
St = S(itest,:) ;
Rt = R(itest,:) ;
S(itest,:) = [] ;
R(itest,:) = [] ;
Tt = length(itest) ;
T = T-Tt ;

save('Retina_stim_resp_network.mat','p','pX','pT','T','Tt','S','R','St','Rt') ;