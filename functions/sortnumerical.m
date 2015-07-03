function files = sortnumerical(infiles)
	indices = [];
	for idx = 1:length(infiles);
		icell = str2num(infiles{idx}(end-5:end-4));
    	if isempty(icell)
    		icell = str2num(infiles{idx}(end-4));
    	end
		indices(idx,1:2) = [icell, idx];
	end
	sorted = sortrows(indices);
	files = infiles(sorted(:,2));
end