load('mabel_reaching_5-4-10')
nU = 45;
spiketimes = {};
isis = {};
unitnames = who('CSPIK*');
dt = 1;
dur = 541187/100;
goodunits = [3 4 6 7 8 9 14 15 16 17 18 19 21 23 27 31 33 34 35 36 37 39 40 41];
for idx = 1:nU
  	%Get data for idxth spike train
  	eval(['spikes = ' unitnames{idx} ';']);
  	spiketimes{idx} = spikes*dt;
  	isis{idx} = diff(spiketimes{idx});
  	below = isis{idx}<500;
  	hist(isis{idx}(below), 100)
  	xlabel('inter-spike interval')
  	ylabel('count')
  	titlestr = ['unit: ' num2str(idx) ', no. spikes: ' num2str(length(isis{idx}))...
  	 ', min(isi): ' num2str(min(isis{idx}))...
  	 ' firing rate: ' num2str(length(isis{idx})/dur)];
  	if ismember(idx, goodunits)
  		titlestr = [titlestr, ', ''good'''];
  	else
  		titlestr = [titlestr, ', ''bad'''];
  	end 
  	title(titlestr)
	saveplot(gcf, ['./results_isis/unit_' num2str(idx) '.eps'])
end


