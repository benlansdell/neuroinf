%Select which units to analyze
goodunits = [4,7,14,15,17,20,24,36,41];

global RefreshRate;
RefreshRate = 100;              %Stimulus refresh rate (Hz)

ds = 0.001;                     %Spike time resolution (seconds)
dt = ds*RefreshRate;            %Time scale spike times are resovled at,
                                % relative to stim timescale
datafile = './mabel_reaching_5-4-10.mat';
nS = 4;                         %no. stim components
frames = 200;                    %no. stim frames 
nF = 2*frames+1;
p = nF*nS;                      %no. stim parameters 
binsize = 1/RefreshRate;
standardize = 0;
%Prepare data
[proc, proc_withheld] = preprocess(datafile, binsize, dt, frames, standardize, goodunits);    

STAs = {};
for ii = 1:length(goodunits)
	icell = goodunits(ii);
    disp(['STA for unit ' num2str(icell)]);
    stim = proc.stim;
    stim = stim/p;
    sptrain = proc.spiketrain(:,ii);
    %Stack stimulus
    stacked = proc.stacked;
    stacked = stacked/p;
    %Compute STA
    STAs{ii} = stacked'*sptrain/sum(sptrain)-mean(stacked,1)'; 
    %Reformat
    STAs{ii} = reshape(STAs{ii},nF,[]);
end