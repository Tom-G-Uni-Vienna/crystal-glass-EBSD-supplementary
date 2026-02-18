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

%% varying minsize
% Processing thresholds (how these were defined is discussed in the paper)
ci_thresh=0.1; % confidence index threshold
fit_thresh=1.4; % fit threshold
grainthresh=15; % grain angle threshold
smooth_runs=15; % N of boundary smooth repetitions

% List of min size thresholds to carry out
minsizethresh=[1 2 3 4 5 6 7 8 9 10];
syze=size(minsizethresh);

% Pre-define empty cell variables of the correct size
minsize_results=cell(size(datastruct));
minsize_value=minsize_results;

minsize_grains=cell(size(minsizethresh));
majoraxM=minsize_grains;
majoraxD=minsize_grains;
minoraxM=minsize_grains;
minoraxD=minsize_grains;

N_density_Mt_true=minsize_grains;
area_Tmt_grains=N_density_Mt_true;
area_Di_grains=N_density_Mt_true;
area_Esk_grains=N_density_Mt_true;
area_notindexed_grains=N_density_Mt_true;
GBlength_Mt=N_density_Mt_true;
GBlength_Di=N_density_Mt_true;

% find size of datastruct (i.e. the input data structure)
size_datastruct=size(datastruct);

% Start looped calculation (looping over each scan and both 
% non-standardised and standardised scans)
for scan_num=1:size_datastruct(1)

fprintf(['\nloading map number ' num2str(scan_num) ' of ' num2str(size_datastruct(1)) '\n'])  

    for proc_type=1:size_datastruct(2) 
    
    fprintf(['\nloading map processing type ' num2str(proc_type) ' of ' num2str(size_datastruct(2)) '\n'])
    
    ebsd=datastruct(scan_num,proc_type).data;
    true_area=datastruct(scan_num,proc_type).true_area;
    
    fprintf('map loaded! cleaning and calculating initial grains...\n') 

    % Clean scans by setting all points outside of quality thresholds to not
    % indexed
    ebsd(ebsd.ci<ci_thresh).phase=0;
    ebsd(ebsd.fit>fit_thresh).phase=0;
    
    % Calc initial grains
    [grains, ebsd.grainId, ebsd.mis2mean]=calcGrains(ebsd,'angle', grainthresh*degree);
    
    ebsd_orig=ebsd;
    grains_orig=grains;
    
            % for a given scan + scan type, calculate grains for all
            % desired min grain size thresholds
            for kk=1:length(minsizethresh)
                
                ebsd=ebsd_orig;
                grains=grains_orig;
    
                fprintf(['Evaluating minsize number ' num2str(kk) ' of ' num2str(length(minsizethresh)) '... \n']) 
                
                % Identify grains with size below min size threshold, and delete the EBSD
                % pixels belonging to them
                toRemove = grains(grains.grainSize<minsizethresh(kk));
                ebsd(toRemove)=[];
                
                clear toRemove grains 
    
                % Calculate final grains after removing pixels belonging to
                % grains that are below the min threshold
                [grains2, ebsd.grainId, ebsd.mis2mean]=calcGrains(ebsd,'angle', grainthresh*degree); 
                
                % Smooth boundaries
                grains2=smooth(grains2,smooth_runs);
    
                % Save result in a cell array for later
                minsize_grains{kk}=grains2;  

                % % % Obtain microstructural parameters % % % % %

                % Tmt number density
                N_density_Mt_true{kk}=length(grains2('Magnetite'))/(true_area/1e6);    
                 
                % Grain areas 
                area_Tmt_grains{kk}=sum(grains2('Magnetite').area);
                area_Di_grains{kk}=sum(grains2('Diopside').area);
                if max(max(ebsd.phase,[],2))>2
                area_Esk_grains{kk}=sum(grains2('Eskolaite').area);
                else
                area_Esk_grains{kk}=0;
                end
                area_notindexed_grains{kk}=sum(grains2('notIndexed').area);

                % % % Grain sizes (from ellipsoids) % % % %
                % MTEX provides the semi-major and semi-minor axes of the ellipse
                % -> these values are multiplied by 2 to obtain total length
                % and width of ellipse, more commonly used as "grain size" in
                % petrology etc.
                
                [~,majoraxM{kk},minoraxM{kk}] = grains2('Magnetite').fitEllipse;
                majoraxM{kk}=majoraxM{kk}.*2;
                minoraxM{kk}=minoraxM{kk}.*2;
                
                [~,majoraxD{kk},minoraxD{kk}] = grains2('Diopside').fitEllipse;
                majoraxD{kk}=majoraxD{kk}.*2; 
                minoraxD{kk}=minoraxD{kk}.*2;
    
                % Boundary lengths
                GBlength_Di{kk}=sum(grains2('Diopside').boundary.segLength);
                GBlength_Mt{kk}=sum(grains2('Magnetite').boundary.segLength);
    
                clear grains2 ebsd
    
            end
    
        clear ebsd_orig grains_orig
        
        % Create a data structure the size of the input minsizethresh list
        % containing the following data, under fieldnames as defined here
        result=struct('minsize_grains',minsize_grains,...
            'majoraxM',majoraxM,'majoraxD',majoraxD,...
            'minoraxM',minoraxM,'minoraxD',minoraxD,...
            'N_density_Mt_true',N_density_Mt_true,...
            'area_Tmt_grains',area_Tmt_grains,'area_Di_grains',area_Di_grains,...
            'area_Esk_grains',area_Esk_grains,'area_notindexed_grains',area_notindexed_grains,...
            'GBlength_Mt',GBlength_Mt,'GBlength_Di',GBlength_Di)  ;    
    
            minsize_results{scan_num,proc_type}=result;
            minsize_value{scan_num,proc_type}=minsizethresh;
            
       clear result
    end

end

fprintf('\nsaving results...\n') 

% Add new fields to the existing datastructure, containing the results of 
% calculating grains for each scan using the list of input minimum 
% grainsize threshold values
[datastruct.minsize_result] = minsize_results{:};
[datastruct.minsize_value] = minsize_value{:};

cd(savelocation)

% Save the updatated data structure as a new variable (large files need the
% option '-v7.3' to allow saving)
save('ETK5a5b_allscans_distort_minsizethresh.mat','datastruct','-v7.3')

fprintf('done!\n\n') 

