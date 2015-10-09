wd = pwd;
if (~isdeployed)
	addpath([pwd '/functions']);
	addpath([pwd '/scripts']);
	addpath_recurse([pwd '/GLM_Algorithm_Functions']);
	addpath_recurse(['~/matlab/schmidt/']);
	addpath_recurse(['~/matlab/chronux/']);
	addpath('~/matlab/plot2svg');
end
