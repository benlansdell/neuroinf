if ismac
        homedir = '/Users/';
else
        homedir = '/home/';
end

addpath([homedir 'lansdell/projects/neuroinf/code']);
addpath([homedir 'lansdell/projects/neuroinf/data']);
addpath([homedir 'lansdell/projects/neuroinf/functions']);
addpath([homedir 'lansdell/projects/neuroinf/eval']);
addpath([homedir 'lansdell/projects/neuroinf/fitting']);
addpath([homedir 'lansdell/projects/neuroinf/preprocess']);
addpath([homedir 'lansdell/projects/neuroinf/models']);
addpath([homedir 'lansdell/projects/neuroinf/scripts']);
addpath_recurse([homedir 'lansdell/projects/neuroinf/code']);


