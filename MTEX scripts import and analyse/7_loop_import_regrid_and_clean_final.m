%% Specify File Names

% Enter the location of the supplementary data folder
supplementary_location =...
    'Z:\Crystal clustering project RESEARCH\Papers\EBSD vs BSE paper';

% Path to files (automatically derived as long as supp. mat. structure
% remains unchanged)
pname = [ supplementary_location...
    '\FINAL Supplementary\FINAL scripts & data\Raw data'];

% Output path (automatically derived as long as supp. mat. structure
% remains unchanged)
savelocation = [ supplementary_location...
    '\FINAL Supplementary\FINAL scripts & data\MATLAB data'];

% Location of the main import & processing script (automatically derived 
% as long as supp. mat. structure remains unchanged)
scriptlocation = [ supplementary_location...
    '\FINAL Supplementary\FINAL scripts & data\MTEX scripts import and analyse'];

% List of the files to be imported
fileNameList={
    '\ETK_5a_scan01_Rescan3'
    '\ETK_5a_scan05_Rescan3'
    '\ETK_5b_scan06_Rescan3'
    '\ETK_5b_scan07_Rescan'
    };

% List of the types of scan - unprocessed ('') and standardised
% ('_CIS_FitS')
preprocessList={
    ''
    '_CIS_FitS'
    };

%% Specify Crystal & Specimen Symmetries

% Crystal symmetry list
CS = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [8.396 8.396 8.396], 'mineral', 'Magnetite', 'color', 'red'),...
  crystalSymmetry('2/m', [9.585 8.776 5.26], [90,106.85,90]*degree, 'X||a*', 'Y||b*', 'Z||c', 'mineral', 'Diopside', 'color', [0.8 0.8 0.8]),...
  crystalSymmetry('-3m', [4.95 4.95 13.58], 'X||a', 'Y||b*', 'Z||c*', 'mineral', 'Eskolaite'),...
  crystalSymmetry('-1', [8.152 12.834 7.079], [93.49,116.13,90.4]*degree,...
  'X||a*', 'Z||c', 'mineral', 'Labradorite', 'color', [0.85 0.65 0.13])...
  };

% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

SS = specimenSymmetry('-1');

%%
% Processing thresholds (how these were defined is discussed in the paper)
ci_thresh=0.1; % confidence index threshold
fit_thresh=1.4; % fit threshold
minsizethresh=5; % minimum possible grain size in pixels
grainthresh=15; % grain angle threshold
smooth_runs=15; % number of smoothing repetitions for boundaries

% Instructions to main script regarding re-gridding:
% Values > 1 = "don't calculate", value of 0 = import but no regrid, all
% other values = step size for re-gridding in µm.
regrid_table= [ 0.2   0.0   1.0    0.0  ];

% Measurements of true vs apparent scan corner shifts from SE image (in µm)
inputshifts=[   316.93	5.41	-316.17	1.47	315.91	 2.21
                304.05	4.04	-305.13	2.10	303.87	 1.78
                326.63	6.31	-329.41	2.56	325.83	 3.49
                510.99	7.25	-504.59	1.02	506.82	-5.71  ];

% % % Run the processing script for the CI threshold chosen (0.1)

fprintf (['\nProcessing with CI threshold = ' num2str(ci_thresh) ': \n'])

cd(scriptlocation)

    run('do_import_regrid_and_clean_final.m')

cd(savelocation)

save('ETK5a5b_allscans_distort_regrid_final.mat','datastruct','regrid_table','ci_thresh','fit_thresh',...
                'minsizethresh', 'grainthresh', 'smooth_runs', '-v7.3')

clear datastruct

% % % Run the processing script for CI threshold 0.075

ci_thresh=0.075;

fprintf (['Processing with CI threshold = ' num2str(ci_thresh) ': \n\n'])

cd(scriptlocation)

    run('do_import_regrid_and_clean_final.m')

datastruct075=datastruct;

cd(savelocation)

save('ETK5a5b_allscans_distort_regrid_final_075.mat','datastruct075','regrid_table','ci_thresh','fit_thresh',...
                'minsizethresh', 'grainthresh', 'smooth_runs', '-v7.3')

clear  datastruct datastruct075

% % % Run the processing script for CI threshold 0.125

ci_thresh=0.125;

fprintf (['Processing with CI threshold = ' num2str(ci_thresh) ': \n\n'])

cd(scriptlocation)

    run('do_import_regrid_and_clean_final.m')

datastruct125=datastruct;

fprintf ('Saving all data... \n')

cd(savelocation)

save('ETK5a5b_allscans_distort_regrid_final_125.mat','datastruct125','regrid_table','ci_thresh','fit_thresh',...
                'minsizethresh', 'grainthresh', 'smooth_runs', '-v7.3')

clear  datastruct  datastruct125

fprintf ('Done! \n\n')