%-Configfile for preprocessfmri.m
%- for spm 12
%__________________________________________________________________________

paralist.spmversion = 'spm12'; 

%-Please specify parallel or nonparallel
%-e.g. for preprocessing and individualstats, set to 1 (parallel)
%-for groupstats, set to 0 (nonparallel)
paralist.parallel = '1';

%-Subject list
paralist.subjectlist = '';
%---- example ----
%- PID, visit, session
%- 7014, 1 ,1 

%-Run list
paralist.runlist = '';
%----- example -------
% task-run01
% task-run02

%-The entire preprocessing to be completed
%-Choose from: 'swar',  'swavr', 'swaor', 'swgcar',  'swgcavr', 'swgcaor'
%-             'swfar', 'swfavr', 'swfaor', 'swgcfar', 'swgcfavr', 'swgcfaor'
%-"s" is smoothing; "w" is normalization; "a": slice timing correction ; "r": realignment
%-"c" is coregistration; "g" is use segmented t1 images while
%coregistration
%-"v" is the 1st version and "o" is the 2nd version of VolRepair pipeline
%-"f" is for fmri images that were acquired flipped
paralist.pipeline = '';

% I/O parameters
% - Raw data directory
paralist.rawdatadir = '';

% - Project directory - output of the preprocessing will be saved in the
% data/imaging folder of the project directory
paralist.projectdir = '';

% - Output directory name
paralist.outputdirname = '';

% fMRI parameters
% - spm8 batch templates location
paralist.batchtemplatepath = '';

% - prefix for the unnormalized file name. for scsnl lab. this value is
% usually empty
paralist.inputimgprefix = '';

% - TR value
paralist.trval = ;

% - Custom slice timing
paralist.customslicetiming = 1;

% - file that specifies the slice timing, in case customslicetiming is
% equal to 1
paralist.slicetimingfile = '';

% - smoothing kernel
paralist.smoothwidth       = [6 6 6];

% - bounding box
paralist.boundingboxdim     = [-90 -126 -72; 90 90 108];

% Coregistration Parameters:
%-Additinal subject list for swgc** pipelines due to better SPGR quality,
% one-to-one matched to paralist.SubjectList
paralist.spgrsubjectlist = ' ';
%%---- example of spgrsubjectlist.csv ------
% PID, visit, session
% 7014, 1 ,1 


% name of the skullstriped T1w volume:
% 'spgr' for original image
% 'skullstrip_spgr_spm12' for brains stripped with SPM12 pipeline
% 'watershed_spgr' for brains stripped with mri_watershed

paralist.spgrfilename = 'spgr';
