%% Specify File Names and load uncleaned data

% Enter the location of the supplementary data folder
supplementary_location =...
    'Z:\Crystal clustering project RESEARCH\Papers\EBSD vs BSE paper';

% Path to files (automatically derived as long as supp. mat. structure
% remains unchanged)
loadlocation = [ supplementary_location...
    '\FINAL Supplementary\FINAL scripts & data\MATLAB data'];

% Output path (automatically derived as long as supp. mat. structure
% remains unchanged)
savelocation = loadlocation;

fprintf('\nLoading MATLAB variable...\n') 

cd(loadlocation)
load('ETK5a5b_allscans_distort.mat')

fprintf('MATLAB variable loaded!\n\n') 

%% Carry out processing
% Processing thresholds (how these were defined is discussed in the paper)
ci_thresh=0.1; % confidence index threshold
fit_thresh=1.4; % fit threshold
minsizethresh=5; % minimum possible grain size in pixels
grainthresh=15; % grain angle threshold

% Pre-define empty cell variables of the correct size 
ebsdcell=cell(size(datastruct));
grainscell=ebsdcell;

cd(savelocation)

% Here we just loop over all 8 entries in "datastruct" - the ordering of the
% results in the output remains the same as the input (i.e. 1st column
% non-standardised, second column standardised, unchanged scan order).
for ii=1:numel(datastruct)

fprintf(['\ncleaning map number ' num2str(ii) ' of ' num2str(numel(datastruct)) '\n\n']) 

ebsd=datastruct(ii).data;

% Clean scans by setting all points outside of quality thresholds to not
% indexed
ebsd(ebsd.ci<ci_thresh).phase=0;
ebsd(ebsd.fit>fit_thresh).phase=0;

fprintf('doing calcgrains one...\n') 

[grains, ebsd.grainId, ebsd.mis2mean]=calcGrains(ebsd,'angle', grainthresh*degree);

% Identify grains with size below min size threshold, and delete the EBSD
% pixels belonging to them
toRemove = grains(grains.grainSize<minsizethresh);
ebsd(toRemove)=[];

clear grains toRemove 

fprintf('doing calcgrains two...\n') 

[grains2, ebsd.grainId, ebsd.mis2mean]=calcGrains(ebsd,'angle', grainthresh*degree);

ebsdcell{ii}=ebsd;
grainscell{ii}=grains2;

clear ebsd grains2

end

% Add new fields to the existing datastructure, containing the cleaned EBSD
% and calculated grains data
[datastruct.cleanebsd] = ebsdcell{:};
[datastruct.grains2] = grainscell{:};

fprintf('\nsaving final result...\n') 

% Save the updatated data structure as a new variable (large files need the
% option '-v7.3' to allow saving)
save('ETK5a5b_allscans_distort_clean.mat','datastruct','-v7.3')

fprintf('done!\n') 