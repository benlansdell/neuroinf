load('mabel_reaching_5-4-10')
nU = 45;
spiketimes = {};
isis = {};
unitnames = who('CSPIK*');
dt = 1;
dur = 541187/100;
goodunits = [3 4 6 7 8 9 14 15 16 17 18 19 21 23 27 31 33 34 35 36 37 39 40 41];
threshold = 100;
rp = 5; %In milliseconds
for idx = 1:nU
  	%Get data for idxth spike train
  	eval(['spikes = ' unitnames{idx} ';']);
  	spiketimes{idx} = spikes*dt;
  	isis{idx} = diff(spiketimes{idx});
  	below = isis{idx}<threshold;
    rpv = sum(isis{idx}<rp);
  	%hist(log(isis{idx}(below)), 100)
    c_min = -1.2; c_max = 2.2;
    edges = 10.^(c_min:0.1:c_max); 
    h = histc(isis{idx}(below), edges);
    centers = sqrt(edges(1:end-1).*edges(2:end));
    bar(h)
    %# fix the x-labels, x-axis extents
    xlim([0.5,length(centers)+0.5])
    tickpts = [.1, 1, 10, 100];
    dcenters = log10(tickpts);
    lcenters = log10(centers);
    ticks = .5+length(centers)*(dcenters-lcenters(1))/(lcenters(end)-lcenters(1));
    set(gca,'xtick',ticks)
    strngs = cellfun(@(x) ['10^{' num2str(log10(x)) '}'], num2cell(tickpts), 'UniformOutput', false);    
    set(gca,'xticklabel',strngs);

  	xlabel('inter-spike interval (ms)')
  	ylabel('count')
  	titlestr = ['unit: ' num2str(idx) ', no. spikes: ' num2str(length(isis{idx}))...
  	 ', min(isi): ' num2str(min(isis{idx}))...
  	 ' firing rate: ' num2str(length(isis{idx})/dur)];
  	if ismember(idx, goodunits)
  		titlestr = [titlestr, ', ''good'''];
  	else
  		titlestr = [titlestr, ', ''bad'''];
  	end 

    %Compute refractory period violations
    N = length(isis{idx});
    T = dur;
    RP = 0.005;
    RPV = rpv;
    [ev,lb,ub] = rpv_contamination(N, T, RP, RPV );

    titlestr = [titlestr, ', contam rate: ' num2str(ev)]
  	title(titlestr)

  	saveplot(gcf, ['./results_isis/unit_' num2str(idx) '_log.eps'])
end









