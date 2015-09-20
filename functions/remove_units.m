function [processed, processed_withheld] = remove_bad_units(indices, proc, proc_withheld)
	processed = proc;
	processed.spikes = processed.spikes(indices);
	processed.spiketrain = processed.spiketrain(:,indices);
	processed.unitnames = processed.unitnames(indices);

	processed_withheld = proc_withheld;
	processed_withheld.spikes = processed_withheld.spikes(indices);
	processed_withheld.spiketrain = processed_withheld.spiketrain(:,indices);
	processed_withheld.unitnames = processed_withheld.unitnames(indices);
end
