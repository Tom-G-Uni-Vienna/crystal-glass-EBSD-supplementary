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

%% Specify Crystal Di Specimen Symmetries

% crystal symmetry
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

% Matrix of instructions regarding re-gridding:
% values > 1 = "don't calculate", value of 0 = import but no regrid, all
% other values = step size for re-gridding. The first map to actually be
% imported has to be without regridding, otherwise the true_area will not
% be correctly calculated (~ line 261)
regrid_table= [ 0.0   9.9   9.9    9.9
                0.2   0.0   0.0    9.9
                0.4   0.4   0.4    9.9
                0.6   0.6   0.6    9.9
                0.8   0.8   0.8    9.9
                1.0   1.0   1.0    0.0  ];

% Measurements of true vs apparent scan corner shifts from SE image (in µm)
inputshifts=[   316.93	5.41	-316.17	1.47	315.91	 2.21
                304.05	4.04	-305.13	2.10	303.87	 1.78                
                326.63	6.31	-329.41	2.56	325.83	 3.49 
                510.99	7.25	-504.59	1.02	506.82	-5.71 ];

% Pre-define empty cell variables of the correct size
% These are of size scan_num * proc_type
data=cell(length(fileNameList),length(preprocessList));
stepsize=data;
mapname=data;
true_area=data;
xspan=data;
yspan=data;

syze=size(regrid_table(:,1));
% These are the same size as the regrid_table
regrid_grains=num2cell(ones(syze));
regrid_ebsd=regrid_grains;
majoraxM=regrid_grains;
majoraxD=regrid_grains;
minoraxM=regrid_grains;
minoraxD=regrid_grains;
smooth_num=regrid_grains;
true_area_regrid=regrid_grains;
stepsize_regrid=regrid_grains;

N_density_Mt_true=regrid_grains;
area_Tmt_grains=N_density_Mt_true;
area_Di_grains=N_density_Mt_true;
area_Esk_grains=N_density_Mt_true;
area_notindexed_grains=N_density_Mt_true;
GBlength_Mt=N_density_Mt_true;
GBlength_Di=N_density_Mt_true;

% Start looped calculation (looping over each scan and both 
% non-standardised and standardised scans)
for scan_num=1:length(fileNameList)  

    fprintf(['\nloading map number ' num2str(scan_num) ' of ' num2str(length(fileNameList)) '\n'])

    for proc_type=1:length(preprocessList) 

        fname = [pname fileNameList{scan_num} preprocessList{proc_type} '.ang'];

        mapname{scan_num,proc_type}=[fileNameList{scan_num} preprocessList{proc_type}];
        
        fprintf(['\nloading map processing type ' num2str(proc_type) ' of ' num2str(length(preprocessList)) '\n'])

        % Generic EBSD import is specifically necessary in MTEX 5.11.2, 
        % along with an updated "loadHelper" script, to avoid the bugged
        % version of loadEBSD_ang of that version. See Supplementary Text.
        ebsd_orig = EBSD.load(fname,CS,SS,'interface','generic' ...
              , 'radians', 'ColumnNames',...
              { 'Euler 1' 'Euler 2' 'Euler 3' 'x' 'y' 'iq' 'ci'  'Phase' 'fit'},...
              'Columns', [1 2 3 4 5 6 7 8 10], 'Bunge');
    
        % Uni Wien reference frame setup correction
        rot1 = rotation.byAxisAngle(zvector,90*degree);
        ebsd_orig = rotate(ebsd_orig,rot1,'keepXY');
        rot2 = rotation.byAxisAngle(xvector,180*degree);
        ebsd_orig = rotate(ebsd_orig,rot2,'keepEuler');

        % Obtain undistorted size
        ext=ebsd_orig.extent;
        xspan{scan_num,proc_type}=(ext(2)-ext(1));
        yspan{scan_num,proc_type}=(ext(4)-ext(3));

        stepsize{scan_num,proc_type}=2*max(ebsd_orig.unitCell(:,1));
    
        fprintf 'map loaded! starting cleaning...\n'         

        % Now regrid the map according to the regrid table instructions
        % (same table applies regardless of which proc_type value)
        for pp=1:length(regrid_table(:,scan_num))

            fprintf(['check regrid instruction ' num2str(pp) ' of ' num2str(length(regrid_table(:,scan_num))) '\n'])

            if regrid_table(pp,scan_num)>1
                % Return empty output here for the final structure
                data{scan_num,proc_type}=struct;

                fprintf 'skipping this one...\n'                
                continue
            else
            end    
            
            ebsd=ebsd_orig;

            % % % % % Do regridding % % % % %
            fprintf (['initial size ' num2str(2*max(ebsd.unitCell(:,1))) 'µm\n'])
            
            if regrid_table(pp,scan_num)==0
            fprintf 'calculating without regrid! \n'   
            else
            fprintf (['regridding to ' num2str(regrid_table(pp,scan_num)) ' um! \n'])

            % Calculate new unitCell size   
            scaling=regrid_table(pp,scan_num)/(2*max(ebsd.unitCell(:,1)));
            newcell=ebsd.unitCell.*scaling;     
              
            % Define the EBSD data set on this new grid
            ebsd = gridify(ebsd,'unitCell',newcell,'nearest');

            end 

            % Make sure that the EBSD variable is not "gridified"
            if length(max(ebsd.x)) > 1
            ebsd=ebsd(:);
            else
            end

            stepsize_regrid{pp}=2*max(ebsd.unitCell(:,1));

            % % % % Start of distort to true map shape based on SE imaging

            % First, centre the EBSD x-y data (middle = 0,0)
            ebsd.x=ebsd.x-((max(ebsd.x)+min(ebsd.x))/2);
            ebsd.y=ebsd.y-((max(ebsd.y)+min(ebsd.y))/2);

            % Determine the ID number of the EBSD pixels at each corner of the
            % scan
            topleft_id=findByLocation(ebsd,[-max(ebsd.x),max(ebsd.y)]);
            topright_id=findByLocation(ebsd,[max(ebsd.x),max(ebsd.y)]);
            botright_id=findByLocation(ebsd,[max(ebsd.x),min(ebsd.y)]);
            botleft_id=findByLocation(ebsd,[-max(ebsd.x),min(ebsd.y)]);

            % The order of the corners is always
            % [bottom left, bottom right, top right, top left]
            % in the following
    
            % Obtain the X (du) and Y (dv) shifts of the scan corners with
            % respect to the nominal shape of the scan as defined during scan
            % acquisition, using the inputshifts above, in µm.                    
            du=[inputshifts(scan_num,4), inputshifts(scan_num,5)-(2*max(ebsd.x)) , inputshifts(scan_num,1)-(2*max(ebsd.x)) , 0];
            dv=[inputshifts(scan_num,3)-(2*-max(ebsd.y)), inputshifts(scan_num,6)+(inputshifts(scan_num,3)-(2*-max(ebsd.y))), inputshifts(scan_num,2), 0];
           
            % normalise the shifts to the nominal scan size
            dunorm=du/max(ebsd.x);
            dvnorm=dv/max(ebsd.y);

            % Define the normalised cordinates of the four corners of the
            % nominal scan (u = X, v = Y)            
            unorm=[-1,1,1,-1];
            vnorm=[-1,-1,1,1];
            
            % Obtain the normalised positions of the scan corners after
            % deformation applied, by summing original positions plus shifts
            unorm=unorm+dunorm;
            vnorm=vnorm+dvnorm;
            
            % Obtain the normalised x-y coordinates of every pixel in the
            % undeformed map
            xnorm=ebsd.x/max(ebsd.x);
            ynorm=ebsd.y/max(ebsd.y);
            
            % Create empty variables which will store the nomralised x and y
            % positions after deformation
            x2=zeros(length(xnorm),1);
            y2=x2;

            % Calculate the normalised x and y positions of each pixel after
            % deformation, based on the positions of the corners
            for q=1:length(xnorm)
                unew=   0.25*((1-xnorm(q))*(1-ynorm(q))*unorm(1))+...
                        0.25*((1+xnorm(q))*(1-ynorm(q))*unorm(2))+...
                        0.25*((1+xnorm(q))*(1+ynorm(q))*unorm(3))+...
                        0.25*((1-xnorm(q))*(1+ynorm(q))*unorm(4));
            
                vnew=   0.25*((1-xnorm(q))*(1-ynorm(q))*vnorm(1))+...
                        0.25*((1+xnorm(q))*(1-ynorm(q))*vnorm(2))+...
                        0.25*((1+xnorm(q))*(1+ynorm(q))*vnorm(3))+...
                        0.25*((1-xnorm(q))*(1+ynorm(q))*vnorm(4));
            
                x2(q)=unew;
                y2(q)=vnew;
            end
            
            % Reverse normalisation of the dataset to obtain new coordinates of
            % pixels in deformed map
            ebsd.x=x2*max(ebsd.x);
            ebsd.y=y2*max(ebsd.y); 

            % % % End of distort to true map shape based on SE imaging % % %           
           
            % Calculate the true area of the distorted scan using the x and y
            % coordinates of the four corners, identified from the id
            % numbers obtained earlier
            eyedees=[topleft_id,topright_id,botright_id,botleft_id];

            if regrid_table(pp,scan_num)==0
            true_area{scan_num,proc_type}=polyarea(ebsd(eyedees).x , ebsd(eyedees).y);
            end

            % % % % Cleaning and grain calculation % % % % %
            fprintf ('beginning cleaning... \n')

            % Clean scans by setting all points outside of quality thresholds to not
            % indexed
            ebsd(ebsd.ci<ci_thresh).phase=0;
            ebsd(ebsd.fit>fit_thresh).phase=0;
            
            fprintf ('calculating grains... \n')

            [grains, ebsd.grainId, ebsd.mis2mean]=calcGrains(ebsd,'angle', grainthresh*degree);
 
            % Identify grains with size below min size threshold, and delete the EBSD
            % pixels belonging to them            
            toRemove = grains(grains.grainSize<minsizethresh);
            ebsd(toRemove)=[];
            
            clear grains toRemove

            fprintf ('calculating grains2... \n')

            [grains2, ebsd.grainId, ebsd.mis2mean]=calcGrains(ebsd,'angle', grainthresh*degree);
            
            % Smooth boundaries
            grains2=smooth(grains2,smooth_runs);

            fprintf ('calculating microstructural parameters... \n')
            % % % Obtain microstructural parameters % % % % %

            % Tmt number density
            N_density_Mt_true{pp}=length(grains2('Magnetite'))/(true_area{scan_num,proc_type}/1e6);

            % Grain areas 
            area_Tmt_grains{pp}=sum(grains2('Magnetite').area);
            area_Di_grains{pp}=sum(grains2('Diopside').area);
            if max(max(ebsd.phase,[],2))>2
            area_Esk_grains{pp}=sum(grains2('Eskolaite').area);
            else
            area_Esk_grains{pp}=0;
            end
            area_notindexed_grains{pp}=sum(grains2('notIndexed').area);

            % % % Grain sizes (from ellipsoids) % % % %
            % MTEX provides the semi-major and semi-minor axes of the ellipse
            % -> these values are multiplied by 2 to obtain total length
            % and width of ellipse, more commonly used as "grain size" in
            % petrology etc.
            [~,majoraxM{pp},minoraxM{pp}] = grains2('Magnetite').fitEllipse;
            majoraxM{pp}=majoraxM{pp}.*2;
            minoraxM{pp}=minoraxM{pp}.*2;

            [~,majoraxD{pp},minoraxD{pp}] = grains2('Diopside').fitEllipse;
            majoraxD{pp}=majoraxD{pp}.*2;     
            minoraxD{pp}=minoraxD{pp}.*2;

            % Boundary lengths
            GBlength_Di{pp}=sum(grains2('Diopside').boundary.segLength);
            GBlength_Mt{pp}=sum(grains2('Magnetite').boundary.segLength);

            regrid_ebsd{pp}=ebsd;

            regrid_grains{pp}=grains2;       
            
        % Create a data structure the size of the input regrid_table
        % containing the following data, under fieldnames as defined here,
        % and add to cell array "data"
            data{scan_num,proc_type}=struct('stepsize_regrid',stepsize_regrid,...    
            'regrid_grains',regrid_grains,'majoraxM',majoraxM,'majoraxD',majoraxD,...
            'minoraxM',minoraxM,'minoraxD',minoraxD,...
            'N_density_Mt_true',N_density_Mt_true,...           
            'area_Tmt_grains',area_Tmt_grains,'area_Di_grains',area_Di_grains,...
            'area_Esk_grains',area_Esk_grains,'area_notindexed_grains',area_notindexed_grains,...
            'GBlength_Mt',GBlength_Mt,'GBlength_Di',GBlength_Di)  ;
    
            clear ebsd grains2 x2 y2 xnorm ynorm dunorm dvnorm unorm vnorm unew vnew topleft_id topright_id botright_id botleft_id eyedees ext

        end      

        clear ebsd_orig 
    end

    % Re-blank these cell arrays of size(regrid_table) before going on to the next scan
    regrid_grains=num2cell(ones(syze));
    regrid_ebsd=regrid_grains;
    majoraxM=regrid_grains;
    majoraxD=regrid_grains;
    smooth_num=regrid_grains;
    true_area_regrid=regrid_grains;
    stepsize_regrid=regrid_grains;
    
    N_density_Mt_true=regrid_grains;    
    area_Tmt_grains=N_density_Mt_true;
    area_Di_grains=N_density_Mt_true;
    area_Esk_grains=N_density_Mt_true;
    area_notindexed_grains=N_density_Mt_true;
    GBlength_Mt=N_density_Mt_true;
    GBlength_Di=N_density_Mt_true;

end

% Create a scan_num x proc_num sized data structure containing the
% following, under fieldnames as defined here. NB data is already itself a
% structure
datastruct=struct('mapname',mapname,'stepsize',stepsize,'xspan',xspan,'yspan',yspan,'true_area',true_area,'data',data);

cd(savelocation)

fprintf ('\nsaving result...\n')

% Save the data structure as a new variable (large files need the
% option '-v7.3' to allow saving)
save('ETK5a5b_allscans_distort_regrid.mat','datastruct','regrid_table','ci_thresh','fit_thresh',...
                'minsizethresh', 'grainthresh', 'smooth_runs', '-v7.3')

fprintf ('done!\n')