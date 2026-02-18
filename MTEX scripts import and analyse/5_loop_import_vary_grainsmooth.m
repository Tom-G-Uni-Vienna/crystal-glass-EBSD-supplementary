%% Specify File Names and load raw data

% Enter the location of the supplementary data folder
supplementary_location =...
    'Z:\Crystal clustering project RESEARCH\Papers\EBSD vs BSE paper';

% Path to saved MATLAB variable (automatically derived as long as supp. mat. structure
% remains unchanged)
loadlocation = [ supplementary_location...
    '\FINAL Supplementary\FINAL scripts & data\MATLAB data'];

% Output path (automatically derived as long as supp. mat. structure
% remains unchanged)
savelocation = loadlocation;

fprintf('\nLoading MATLAB variable...\n') 

cd(loadlocation)
load('ETK5a5b_allscans_distort_clean.mat')

fprintf('MATLAB variable loaded!\n\n') 

%% varying Grainsmooth
% List of smoothing values to carry out
smoothing_values=[1:10,15,20];

% Pre-define empty cell variables of the correct size

% These are of size scan_num * proc_type
smooth_results=cell(size(datastruct));

syze=size(smoothing_values);
% These are the same size as the input list of smoothing values
smooth_grains=cell(syze);
majoraxM=smooth_grains;
majoraxD=smooth_grains;
minoraxM=smooth_grains;
minoraxD=smooth_grains;
smooth_num=smooth_grains;

N_density_Mt_true=smooth_grains;
area_Tmt_grains=N_density_Mt_true;
area_Di_grains=N_density_Mt_true;
area_Esk_grains=N_density_Mt_true;
area_notindexed_grains=N_density_Mt_true;
GBlength_Mt=N_density_Mt_true;
GBlength_Di=N_density_Mt_true;

% Start looped calculation
% Here we just loop over all 8 entries in "datastruct" - the ordering of the
% results in the output remains the same as the input (i.e. 1st column
% non-standardised, second column standardised, unchanged scan order).
for ii=1:numel(datastruct)

fprintf(['\nprocessing grainsmooth for map number ' num2str(ii) ' \n\n']) 

grains_orig=datastruct(ii).grains2;
true_area=datastruct(ii).true_area;

        for kk=1:length(smoothing_values)            
        
            fprintf(['Evaluating smooth number ' num2str(kk) ' of ' num2str(length(smoothing_values)) '. \n']) 

            % Smooth the grain boundaries
            grains_smooth=smooth(grains_orig,smoothing_values(kk));

            smooth_num{kk}=smoothing_values(kk);

            smooth_grains{kk}=grains_smooth;
    
            % % % Obtain microstructural parameters % % % % %

            % Tmt number density
            N_density_Mt_true{kk}=length(grains_smooth('Magnetite'))/(true_area/1e6);
              
            % Grain areas and area fractions 
            area_Tmt_grains{kk}=sum(grains_smooth('Magnetite').area);
            area_Di_grains{kk}=sum(grains_smooth('Diopside').area);
            if max(max(datastruct(ii).grains2.phase,[],2))>2
            area_Esk_grains{kk}=sum(grains_smooth('Eskolaite').area);
            else
            area_Esk_grains{kk}=0;
            end
            area_notindexed_grains{kk}=sum(grains_smooth('notIndexed').area);
            
            % % % Grain sizes (from ellipsoids) % % % %
            % MTEX provides the semi-major and semi-minor axes of the ellipse
            % -> these values are multiplied by 2 to obtain total length
            % and width of ellipse, more commonly used as "grain size" in
            % petrology etc.

            [~,majoraxM{kk},minoraxM{kk}] = grains_smooth('Magnetite').fitEllipse;
            majoraxM{kk}=majoraxM{kk}.*2;
            minoraxM{kk}=minoraxM{kk}.*2;

            [~,majoraxD{kk},minoraxD{kk}] = grains_smooth('Diopside').fitEllipse;
            majoraxD{kk}=majoraxD{kk}.*2;
            minoraxD{kk}=minoraxD{kk}.*2;

            % Boundary lengths
            GBlength_Di{kk}=sum(grains_smooth('Diopside').boundary.segLength);
            GBlength_Mt{kk}=sum(grains_smooth('Magnetite').boundary.segLength);

            clear grains_smooth

        end

    % Create a data structure the size of the input smooth values list
    % containing the following data, under fieldnames as defined here
    result=struct('smooth_grains',smooth_grains,'smooth_num',smooth_num,...
        'majoraxM',majoraxM,'majoraxD',majoraxD,...
        'minoraxM',minoraxM,'minoraxD',minoraxD,...
        'N_density_Mt_true',N_density_Mt_true,...
        'area_Tmt_grains',area_Tmt_grains,'area_Di_grains',area_Di_grains,...
        'area_Esk_grains',area_Esk_grains,'area_notindexed_grains',area_notindexed_grains,...
        'GBlength_Mt',GBlength_Mt,'GBlength_Di',GBlength_Di)  ;    

    smooth_results{ii}=result;

    clear result    

end

fprintf('\nsaving smooth results...\n'); 

% Add new field to the existing datastructure, containing the results of 
% smoothing each scan using the list of input smoothing values
[datastruct.smooth_result] = smooth_results{:};

cd(savelocation)

% Save the updatated data structure as a new variable (large files need the
% option '-v7.3' to allow saving)
save('ETK5a5b_allscans_distort_grainsmooth.mat','datastruct','-v7.3')

fprintf('done!\n\n') 
