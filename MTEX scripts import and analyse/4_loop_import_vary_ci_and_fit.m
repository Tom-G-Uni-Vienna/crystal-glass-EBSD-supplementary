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
load('ETK5a5b_allscans_distort.mat')

fprintf('MATLAB variable loaded!\n\n') 
%% varying Ci and fit
% construct matrix of CI and fit values to calculate maps for
ci_vary=[0.05:0.025:0.15]';
ci_ignore=ones(length(ci_vary),1).*-2;
ci_const=ones(length(ci_vary),1).*0.1;

fit_vary=[1.2:0.1:1.6]';
fit_ignore=ones(length(fit_vary),1).*200;
fit_const=ones(length(fit_vary),1).*1.4;

cithresh=horzcat(ci_vary,ci_ignore,ci_vary,ci_const)
fitthresh=horzcat(fit_ignore,fit_vary,fit_const,fit_vary)

syze=size(cithresh);

% Processing thresholds (how these were defined is discussed in the paper)
minsizethresh=5; % minimum possible grain size in pixels
grainthresh=15; % grain angle threshold
smooth_val=15; % number of smoothing repetitions for boundaries

% Pre-define empty cell variables of the correct size

% These are of size scan_num * proc_type
ci_fit_results=cell(size(datastruct));
ci_thresh=ci_fit_results;
fit_thresh=ci_fit_results;

% These are the same size as the input CI and fit matrices
ci_fit_grains=cell(syze);
majoraxM=ci_fit_grains;
majoraxD=ci_fit_grains;
minoraxM=ci_fit_grains;
minoraxD=ci_fit_grains;
N_density_Mt_true=ci_fit_grains;
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

fprintf(['\nprocessing ci_fit_var for map number ' num2str(ii) ' \n\n']) 

ebsd_orig=datastruct(ii).data;
true_area=datastruct(ii).true_area;

inc=0; % incremented number to track progress

    for col=1:syze(2) % columns of the fit & CI matrices

        for row=1:syze(1) % rows of the fit & CI matrices

            inc = inc+1;

            ebsd=ebsd_orig;
        
            fprintf(['column ' num2str(col) ', row ' num2str(row) ' of ' ...
                num2str(syze(2)) 'x' num2str(syze(1)) ' matrix...' ...
                'number ' num2str(inc) ' of 20. \n'])    

            % Clean scans by setting all points outside of quality thresholds to not
            % indexed
            ebsd(ebsd.ci<cithresh(row,col)).phase=0;
            ebsd(ebsd.fit>fitthresh(row,col)).phase=0;
    
            [grains, ebsd.grainId, ebsd.mis2mean]=calcGrains(ebsd,'angle', grainthresh*degree);

            % Identify grains with size below min size threshold, and delete the EBSD
            % pixels belonging to them
            toRemove = grains(grains.grainSize<minsizethresh);
            ebsd(toRemove)=[];
            
            clear grains toRemove

            [grains2, ebsd.grainId, ebsd.mis2mean]=calcGrains(ebsd,'angle', grainthresh*degree); 
            
            % Smooth grain boundaries
            grains2=smooth(grains2,smooth_val);

            % Keep grains result for each CI+Fit pair for later
            ci_fit_grains{row,col}=grains2;

            % % % Obtain microstructural parameters % % % % %

            % Tmt number density
            N_density_Mt_true{row,col}=length(grains2('Magnetite'))/(true_area/1e6);
                      
            % Grain areas and area fractions 
            area_Tmt_grains{row,col}=sum(grains2('Magnetite').area);
            area_Di_grains{row,col}=sum(grains2('Diopside').area);
            if max(max(ebsd.phase,[],2))>2
            area_Esk_grains{row,col}=sum(grains2('Eskolaite').area);
            else
            area_Esk_grains{row,col}=0;
            end
            area_notindexed_grains{row,col}=sum(grains2('notIndexed').area);

            % % % Grain sizes (from ellipsoids) % % % %
            % MTEX provides the semi-major and semi-minor axes of the ellipse
            % -> these values are multiplied by 2 to obtain total length
            % and width of ellipse, more commonly used as "grain size" in
            % petrology etc.
            
            [~,majoraxM{row,col},minoraxM{row,col}] = grains2('Magnetite').fitEllipse;
            majoraxM{row,col}=majoraxM{row,col}.*2;
            minoraxM{row,col}=minoraxM{row,col}.*2;
            
            [~,majoraxD{row,col},minoraxD{row,col}] = grains2('Diopside').fitEllipse;
            majoraxD{row,col}=majoraxD{row,col}.*2;          
            minoraxD{row,col}=minoraxD{row,col}.*2;

            % Boundary lengths
            GBlength_Di{row,col}=sum(grains2('Diopside').boundary.segLength);
            GBlength_Mt{row,col}=sum(grains2('Magnetite').boundary.segLength);

            clear ebsd grains2 % For file size reasons, ebsd maps themselves are not saved

        end    

    end  

    clear ebsd_orig  

    % Create a data structure the size of the input CI+Fit matrices
    % containing the following data, under fieldnames as defined here
    result=struct('ci_fit_grains',ci_fit_grains,...
        'majoraxM',majoraxM,'majoraxD',majoraxD,...
        'minoraxM',minoraxM,'minoraxD',minoraxD,...
        'N_density_Mt_true',N_density_Mt_true,...
        'area_Tmt_grains',area_Tmt_grains,'area_Di_grains',area_Di_grains,...
        'area_Esk_grains',area_Esk_grains,'area_notindexed_grains',area_notindexed_grains,...
        'GBlength_Mt',GBlength_Mt,'GBlength_Di',GBlength_Di)  ;
             
        ci_fit_results{ii}=result;
        ci_thresh{ii}=cithresh;
        fit_thresh{ii}=fitthresh;
        
        clear result

end

fprintf('\nsaving results...\n'); 

% Add new fields to the existing datastructure, containing the row*col
% sized results of processing each scan using the different CI+Fit pairs
% defined by the input CI and fit matrices
[datastruct.cifit_result] = ci_fit_results{:};
[datastruct.cithresh] = ci_thresh{:};
[datastruct.fitthresh] = fit_thresh{:};

cd(savelocation)

% Save the updatated data structure as a new variable (large files need the
% option '-v7.3' to allow saving)
save('ETK5a5b_allscans_distort_cithresh.mat','datastruct','-v7.3')

fprintf('done!\n\n') 