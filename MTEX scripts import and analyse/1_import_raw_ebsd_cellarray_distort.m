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

% Plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

SS = specimenSymmetry('-1');

%% Carry out import

% Measurements of true scan corner locations from SE image (in µm)
% for details of how these are defined, see Supplementary Figures & Text
inputshifts=[   316.93	5.41	-316.17	1.47	315.91	 2.21
                304.05	4.04	-305.13	2.10	303.87	 1.78
                326.63	6.31	-329.41	2.56	325.83	 3.49
                510.99	7.25	-504.59	1.02	506.82	-5.71  ];

% Pre-define empty cell variables of the correct size 
data=cell(length(fileNameList),length(preprocessList));
true_area=data;
stepsize=data;
mapname=data;
xspan=data;
yspan=data;

% Begin the looped import
for scan_num=1:length(fileNameList)  

    fprintf(['\nloading map number ' num2str(scan_num) ' of ' num2str(length(fileNameList)) '\n\n'])

    for proc_type=1:length(preprocessList) 

        fname = [pname fileNameList{scan_num} preprocessList{proc_type} '.ang'];

        mapname{scan_num,proc_type}=[fileNameList{scan_num} preprocessList{proc_type}];
        
        fprintf(['loading map processing type ' num2str(proc_type) ' of ' num2str(length(preprocessList)) '\n\n'])
        
        % Generic EBSD import is specifically necessary in MTEX 5.11.2, 
        % along with an updated "loadHelper" script, to avoid the bugged
        % version of loadEBSD_ang of that version. See Supplementary Text.
        ebsd = EBSD.load(fname,CS,SS,'interface','generic' ...
          , 'radians', 'ColumnNames',...
          { 'Euler 1' 'Euler 2' 'Euler 3' 'x' 'y' 'iq' 'ci'  'Phase' 'fit'},...
          'Columns', [1 2 3 4 5 6 7 8 10], 'Bunge');

        % Uni Wien reference frame setup correction
        rot1 = rotation.byAxisAngle(zvector,90*degree);
        ebsd = rotate(ebsd,rot1,'keepXY');
        rot2 = rotation.byAxisAngle(xvector,180*degree);
        ebsd = rotate(ebsd,rot2,'keepEuler');    

        % Start of distort to true map shape based on SE imaging,
        % this routine is based on the notes at 
        % https://www.continuummechanics.org/finiteelementmapping.html

        % First, centre the EBSD x-y data (middle = 0,0)
        ebsd.x=ebsd.x-((max(ebsd.x)+min(ebsd.x))/2);
        ebsd.y=ebsd.y-((max(ebsd.y)+min(ebsd.y))/2);

        % Copy EBSD variable to keep the original unchanged after centering
        % (needed to calculate nominal size and identify which are the
        % corner pixels later on)
        ebsd_orig=ebsd;
      
        % The order of the corners is always
        % [bottom left, bottom right, top right, top left]
        % in the following

        % Obtain the X (du) and Y (dv) shifts of the scan corners with
        % respect to the nominal shape of the scan as defined during scan
        % acquisition, using the inputshifts above, in µm.
        du=[inputshifts(scan_num,4), inputshifts(scan_num,5)-(2*max(ebsd.x)) , inputshifts(scan_num,1)-(2*max(ebsd.x)) , 0];
        dv=[inputshifts(scan_num,3)-(2*-max(ebsd.y)), inputshifts(scan_num,6)+(inputshifts(scan_num,3)-(2*-max(ebsd.y))), inputshifts(scan_num,2) , 0];
        
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

        % End of distort to true map shape based on SE imaging

        data{scan_num,proc_type}=ebsd;

        % Obtain nominal scan size from the undistorted ebsd_orig variable
        ext=ebsd_orig.extent;
        xspan{scan_num,proc_type}=(ext(2)-ext(1));
        yspan{scan_num,proc_type}=(ext(4)-ext(3));
  
        % Determine the ID number of the EBSD pixels at each corner of the
        % scan
        topleft_id=findByLocation(ebsd_orig,[-max(ebsd_orig.x),max(ebsd_orig.y)]);
        topright_id=findByLocation(ebsd_orig,[max(ebsd_orig.x),max(ebsd_orig.y)]);
        botright_id=findByLocation(ebsd_orig,[max(ebsd_orig.x),min(ebsd_orig.y)]);
        botleft_id=findByLocation(ebsd_orig,[-max(ebsd_orig.x),min(ebsd_orig.y)]);
       
        eyedees=[topleft_id,topright_id,botright_id,botleft_id];
    
        % Calculate the true area of the distorted scan using the x and y
        % coordinates of the four corners
        true_area{scan_num,proc_type}=polyarea(ebsd(eyedees).x , ebsd(eyedees).y);

        % Obtain stepsize variable (warning - specific
        stepsize{scan_num,proc_type}=2*max(ebsd_orig.unitCell(:,1));

        clear ebsd x2 y2 xnorm ynorm dunorm dvnorm unorm vnorm unew vnew ebsd_orig topleft_id topright_id botright_id botleft_id eyedees ext

    end
end

% Create a scan_num x proc_num sized data structure containing the
% following data, under fieldnames as defined here
datastruct=struct('mapname',mapname,'stepsize',stepsize,'xspan',xspan,'yspan',yspan,'true_area',true_area,'data',data);

cd(savelocation)

fprintf('\nsaving result... \n\n')

save('ETK5a5b_allscans_distort.mat','datastruct')

fprintf('done! \n\n')

%%


