%% Specify Plotting Convention and Load data

% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

% Enter the location of the supplementary data folder
supplementary_location =...
    'Z:\Crystal clustering project RESEARCH\Papers\EBSD vs BSE paper';

% Path to load MATLAB variables (automatically derived as long as supp. mat. structure
% remains unchanged)
loadlocation = [ supplementary_location...
    '\FINAL Supplementary\FINAL scripts & data\MATLAB data'];

% Output path (automatically derived as long as supp. mat. structure
% remains unchanged)
savelocation = [ supplementary_location...
    '\FINAL Supplementary\FINAL scripts & data\Figures & tables' ];

% File name for Excel table of microstructural parameters
filenameXLS = 'Pontesilli2019_compareTable.xlsx';

cd(loadlocation)

fprintf('loading MATLAB variables...\n\n') 

% MATLAB prefers to load data to structures and extract variables from them
% instead of loading the variables directly 
% (https://de.mathworks.com/matlabcentral/answers/28676-why-use-x-load-myfile-mat)
input1=load('ETK5a5b_allscans_distort_regrid_final.mat','datastruct');
input2=load('ETK5a5b_allscans_distort_combined_result.mat','resultstruct');

datastruct=input1.datastruct;
resultstruct=input2.resultstruct;

clear input1 input2

fprintf('MATLAB variables loaded!\n\n') 
%% Pontesilli et al. Data

TPManhy =  [ [0.0077 ; 5.31   ; 1.48e-7 ; 0.71e-7 ; 2.56  ; 1174  ], ...
             [0.0068 ; 0.78   ; 0.22e-7 ; 0       ; 0     ; 0     ]];

TPMhy =  [  [0.018  ; 7.20  ; 2.00e-7 ; 1.13e-7 ; 4.08  ; 1083  ], ...
            [0.0157 ; 0.56  ; 0.16e-7 ; 0       ; 0     ; 0     ]];

TPDanhy =  [ [0.5353  ; 53.66 ; 14.91e-7 ; 513.84  ], ...
             [0.1241  ; 5.39  ; 1.50e-7  ; 87.93   ]];

TPDhy =  [ [0.3846  ; 63.24 ; 18.02e-7 ; 398.81  ], ...
           [0.0609  ; 6.41  ; 2.49e-7  ; 92.31   ]];

%% Output Excel table containing microstructural parameter information
% both for combined scans and individual scans

cd(savelocation)

% Helper variable for the combination of data - needed as 
% sample 1 = maps 1+2, sample 2 = maps 3+4 ->
% so the two maps for a given sample in the input variable datastruct are
% at row positions q(sample_num) and q(sample_num)+1.
q=[1,3]; 

% % % % % % % % % Titanomagnetite table % % % % % % % % % % %

% sample
columnletter= {'E';'O'};

% process
startingrow = [ 5 ; 12 ];

writematrix( TPManhy ,filenameXLS,'Sheet','Sheet1','Range',...
    'C5','UseExcel',0);
writematrix( TPManhy ,filenameXLS,'Sheet','Sheet1','Range',...
    'C12','UseExcel',0);

writematrix( TPMhy ,filenameXLS,'Sheet','Sheet1','Range',...
    'M5','UseExcel',0);
writematrix( TPMhy ,filenameXLS,'Sheet','Sheet1','Range',...
    'M12','UseExcel',0);

writematrix( TPDanhy ,filenameXLS,'Sheet','Sheet1','Range',...
    'C22','UseExcel',0);
writematrix( TPDanhy ,filenameXLS,'Sheet','Sheet1','Range',...
    'C27','UseExcel',0);

writematrix( TPDhy ,filenameXLS,'Sheet','Sheet1','Range',...
    'M22','UseExcel',0);
writematrix( TPDhy ,filenameXLS,'Sheet','Sheet1','Range',...
    'M27','UseExcel',0);

for sample_num=1:2
    for proc_type=1:2

    % Temp variable containing combined data for the given sample and
    % processing type (for all three tested CI thresholds)      
    r=resultstruct(sample_num,proc_type);
    
    % Temp variables containing separate data for each scan for the given
    % sample and processing type 
    d1=datastruct(q(sample_num),proc_type);
    d2=datastruct(q(sample_num)+1,proc_type);
    
T = [  [r.afrac_Tmt(1)  ;  r.mean_L10pont_Area_M(1)   ; r.meanGmax_M(1) ; r.Gbatch_Tmt(1) ; r.d_Tmt(1)  ; r.NAtmt(1)  ], ...
       [0               ;  r.std_L10pont_Area_M(1)    ; r.stdGmax_M(1)  ; 0               ; 0           ; 0           ],...
       [r.afrac_Tmt(2)  ;  r.mean_L10pont_Area_M(2)   ; r.meanGmax_M(2) ; r.Gbatch_Tmt(2) ; r.d_Tmt(2)  ; r.NAtmt(2)  ], ...
       [r.afrac_Tmt(3)  ;  r.mean_L10pont_Area_M(3)   ; r.meanGmax_M(3) ; r.Gbatch_Tmt(3) ; r.d_Tmt(3)  ; r.NAtmt(3)  ], ...
       [d1.afrac_Tmt    ;  d1.mean_Lpont10_M          ; d1.mean_Gmax_M  ; d1.Gbatch_Tmt   ; d1.d_Tmt    ; d1.NAtmt    ], ...
       [0               ;  d1.std_Lpont10_M           ; d1.std_Gmax_M   ; 0               ; 0           ; 0           ],...
       [d2.afrac_Tmt    ;  d2.mean_Lpont10_M          ; d2.mean_Gmax_M  ; d2.Gbatch_Tmt   ; d2.d_Tmt    ; d2.NAtmt    ], ...
       [0               ;  d2.std_Lpont10_M           ; d2.std_Gmax_M   ; 0               ; 0           ; 0           ]  ];
    
    writematrix( T ,filenameXLS,'Sheet','Sheet1','Range',...
        [columnletter{sample_num} num2str(startingrow(proc_type))],'UseExcel',0);

    end
end

% % % % % % % % % Clinopyroxene table % % % % % % % % % % %

% process
startingrow = [ 22 ; 27 ];

for sample_num=1:2
    for proc_type=1:2

    % Temp variable containing combined data for the given sample and
    % processing type (for all three tested CI thresholds)     
    r=resultstruct(sample_num,proc_type);
    
    % Temp variables containing separate data for each scan for the given
    % sample and processing type 
    d1=datastruct(q(sample_num),proc_type);
    d2=datastruct(q(sample_num)+1,proc_type);
    
T = [  [r.afrac_Di(1)  ; r.mean_L10pont_Area_D(1)   ; r.meanGmax_D(1) ; r.SVPcpx(1) ], ...
       [0              ; r.std_L10pont_Area_D(1)    ; r.stdGmax_D(1)  ; 0           ],...
       [r.afrac_Di(2)  ; r.mean_L10pont_Area_D(2)   ; r.meanGmax_D(2) ; r.SVPcpx(2) ], ...
       [r.afrac_Di(3)  ; r.mean_L10pont_Area_D(3)   ; r.meanGmax_D(3) ; r.SVPcpx(3) ], ...
       [d1.afrac_Di    ; d1.mean_Lpont10_D          ; d1.mean_Gmax_D  ; d1.SVPcpx   ], ...
       [0              ; d1.std_Lpont10_D           ; d1.std_Gmax_D   ; 0           ],...
       [d2.afrac_Di    ; d2.mean_Lpont10_D          ; d2.mean_Gmax_D  ; d2.SVPcpx   ], ...
       [0              ; d2.std_Lpont10_D           ; d2.std_Gmax_D   ; 0           ]  ];

    
    writematrix( T ,filenameXLS,'Sheet','Sheet1','Range',...
        [columnletter{sample_num} num2str(startingrow(proc_type))],'UseExcel',0);

    end
end