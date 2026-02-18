% Pre-define empty cell variables of the correct size 
data=cell(length(fileNameList),length(preprocessList));
stepsize=data;
mapname=data;
true_area=data;
xspan=data;
yspan=data;

grains2=data;
ebsd=data;

majoraxM=data;
minoraxM=data;
majax10_M=data;
minax10_M=data;
std_Lpont10_M=data;
mean_Lpont10_M=data;
std_majax10_M=data;
mean_majax10_M=data;
std_minax10_M=data;
mean_minax10_M=data;
Gmax_M=data;
std_Gmax_M=data;
mean_Gmax_M=data;

majoraxD=data;
minoraxD=data;
majax10_D=data;
minax10_D=data;
std_Lpont10_D=data;
mean_Lpont10_D=data;
std_majax10_D=data;
mean_majax10_D=data;
std_minax10_D=data;
mean_minax10_D=data;
Gmax_D=data;
std_Gmax_D=data;
mean_Gmax_D=data;

smooth_num=data;
true_area_regrid=data;
stepsize_regrid=data;

N_density_Mt_true=data;
N_Mt_grains=data;

d_Tmt=data;
Gbatch_Tmt=data;

area_Tmt_grains=data;
area_Di_grains=data;
area_Esk_grains=data;
area_notindexed_grains=data;

afrac_Tmt=data;
afrac_Di=data;
afrac_Esk=data;
afrac_notindexed=data;

GBlength_Mt=data;
GBlength_Di=data;

SVPcpx=data;

% Begin the looped import
for scan_num=1:length(fileNameList)  

    fprintf(['\nloading map number ' num2str(scan_num) ' of ' num2str(length(fileNameList)) '\n\n'])

    for proc_type=1:length(preprocessList) 

        fname = [pname fileNameList{scan_num} preprocessList{proc_type} '.ang'];

        mapname{scan_num,proc_type}=[fileNameList{scan_num} preprocessList{proc_type}];
        
        fprintf(['loading map processing type ' num2str(proc_type) ' of ' num2str(length(preprocessList)) '\n'])

        % Generic EBSD import is specifically necessary in MTEX 5.11.2, 
        % along with an updated "loadHelper" script, to avoid the bugged
        % version of loadEBSD_ang of that version. See Supplementary Text.
        ebsd{scan_num,proc_type} = EBSD.load(fname,CS,SS,'interface','generic' ...
              , 'radians', 'ColumnNames',...
              { 'Euler 1' 'Euler 2' 'Euler 3' 'x' 'y' 'iq' 'ci'  'Phase' 'fit'},...
              'Columns', [1 2 3 4 5 6 7 8 10], 'Bunge');
    
        % Uni Wien reference frame setup correction
        rot1 = rotation.byAxisAngle(zvector,90*degree);
        ebsd{scan_num,proc_type} = rotate(ebsd{scan_num,proc_type},rot1,'keepXY');
        rot2 = rotation.byAxisAngle(xvector,180*degree);
        ebsd{scan_num,proc_type} = rotate(ebsd{scan_num,proc_type},rot2,'keepEuler');

        ext=ebsd{scan_num,proc_type}.extent;
        xspan{scan_num,proc_type}=(ext(2)-ext(1));
        yspan{scan_num,proc_type}=(ext(4)-ext(3));

        stepsize{scan_num,proc_type}=2*max(ebsd{scan_num,proc_type}.unitCell(:,1));
    
        fprintf 'map loaded! starting cleaning...\n'         

        % Now regrid the map according to the regrid table instructions
        % (same table applies regardless of which proc_type value)    

            % % % % % Do regridding % % % % %
            fprintf (['initial size ' num2str(2*max(ebsd{scan_num,proc_type}.unitCell(:,1))) 'µm\n'])
            
            if regrid_table(scan_num)==0
            fprintf 'calculating without regrid! \n'   
            else
            fprintf (['regridding to ' num2str(regrid_table(scan_num)) ' um! \n'])
            
            % Calculate new unitCell size
            scaling=regrid_table(scan_num)/(2*max(ebsd{scan_num,proc_type}.unitCell(:,1)));
            newcell=ebsd{scan_num,proc_type}.unitCell.*scaling;     
              
            % Define the EBSD data set on this new grid
            ebsd{scan_num,proc_type} = gridify(ebsd{scan_num,proc_type},'unitCell',newcell,'nearest');

            end  

            % Make sure that the EBSD variable is not "gridified"
            if length(max(ebsd{scan_num,proc_type}.x)) > 1
            ebsd{scan_num,proc_type}=ebsd{scan_num,proc_type}(:);
            end

            stepsize_regrid{scan_num,proc_type}=2*max(ebsd{scan_num,proc_type}.unitCell(:,1));

            % % % % Start of distort to true map shape based on SE imaging
            
            % First, centre the EBSD x-y data (middle = 0,0)
            ebsd{scan_num,proc_type}.x=ebsd{scan_num,proc_type}.x-((max(ebsd{scan_num,proc_type}.x)+min(ebsd{scan_num,proc_type}.x))/2);
            ebsd{scan_num,proc_type}.y=ebsd{scan_num,proc_type}.y-((max(ebsd{scan_num,proc_type}.y)+min(ebsd{scan_num,proc_type}.y))/2);

            % Determine the ID number of the EBSD pixels at each corner of the
            % scan
            topleft_id=findByLocation(ebsd{scan_num,proc_type},[-max(ebsd{scan_num,proc_type}.x),max(ebsd{scan_num,proc_type}.y)]);
            topright_id=findByLocation(ebsd{scan_num,proc_type},[max(ebsd{scan_num,proc_type}.x),max(ebsd{scan_num,proc_type}.y)]);
            botright_id=findByLocation(ebsd{scan_num,proc_type},[max(ebsd{scan_num,proc_type}.x),min(ebsd{scan_num,proc_type}.y)]);
            botleft_id=findByLocation(ebsd{scan_num,proc_type},[-max(ebsd{scan_num,proc_type}.x),min(ebsd{scan_num,proc_type}.y)]);
            
            % The order of the corners is always
            % [bottom left, bottom right, top right, top left]
            % in the following
    
            % Obtain the X (du) and Y (dv) shifts of the scan corners with
            % respect to the nominal shape of the scan as defined during scan
            % acquisition, using the inputshifts above, in µm.
            du=[inputshifts(scan_num,4), inputshifts(scan_num,5)-(2*max(ebsd{scan_num,proc_type}.x)) ,...
                inputshifts(scan_num,1)-(2*max(ebsd{scan_num,proc_type}.x)) , 0];
            dv=[inputshifts(scan_num,3)-(2*-max(ebsd{scan_num,proc_type}.y)),...
                inputshifts(scan_num,6)+(inputshifts(scan_num,3)-(2*-max(ebsd{scan_num,proc_type}.y))),...
                inputshifts(scan_num,2) , 0];

            % normalise the shifts to the nominal scan size
            dunorm=du/max(ebsd{scan_num,proc_type}.x);
            dvnorm=dv/max(ebsd{scan_num,proc_type}.y);

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
            xnorm=ebsd{scan_num,proc_type}.x/max(ebsd{scan_num,proc_type}.x);
            ynorm=ebsd{scan_num,proc_type}.y/max(ebsd{scan_num,proc_type}.y);

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
            ebsd{scan_num,proc_type}.x=x2*max(ebsd{scan_num,proc_type}.x);
            ebsd{scan_num,proc_type}.y=y2*max(ebsd{scan_num,proc_type}.y); 

            % % % End of distort to true map shape based on SE imaging % % %           
           
            % Calculate the true area of the distorted scan using the x and y
            % coordinates of the four corners, identified from the id
            % numbers obtained earlier
            eyedees=[topleft_id,topright_id,botright_id,botleft_id];
            
            true_area{scan_num,proc_type}=polyarea(ebsd{scan_num,proc_type}(eyedees).x , ebsd{scan_num,proc_type}(eyedees).y);

            % % % % Cleaning and grain calculation % % % % %

            fprintf ('beginning cleaning... \n')

            % Clean scans by setting all points outside of quality thresholds to not
            % indexed
            ebsd{scan_num,proc_type}(ebsd{scan_num,proc_type}.ci<ci_thresh).phase=0;
            ebsd{scan_num,proc_type}(ebsd{scan_num,proc_type}.fit>fit_thresh).phase=0;

            fprintf ('calculating grains... \n')

            [grains, ebsd{scan_num,proc_type}.grainId, ebsd{scan_num,proc_type}.mis2mean]=calcGrains(ebsd{scan_num,proc_type},'angle', grainthresh*degree);

            % Identify grains with size below min size threshold, and delete the EBSD
            % pixels belonging to them
            toRemove = grains(grains.grainSize<minsizethresh);
            ebsd{scan_num,proc_type}(toRemove)=[];
            
            clear grains toRemove

            fprintf ('calculating grains2... \n')

            [grains2{scan_num,proc_type}, ebsd{scan_num,proc_type}.grainId, ebsd{scan_num,proc_type}.mis2mean]=calcGrains(ebsd{scan_num,proc_type},'angle', grainthresh*degree);
   
            % Smooth boundaries
            grains2{scan_num,proc_type}=smooth(grains2{scan_num,proc_type},smooth_runs);

            fprintf ('calculating microstructural parameters... \n')

            % % % Obtain microstructural parameters % % % % %

            % N tmt and number density
            N_Mt_grains{scan_num,proc_type}=length(grains2{scan_num,proc_type}('Magnetite'));
            N_density_Mt_true{scan_num,proc_type}=length(grains2{scan_num,proc_type}('Magnetite'))/(true_area{scan_num,proc_type}/1e6);

            % Grain areas and area fractions 
            area_Tmt_grains{scan_num,proc_type}=sum(grains2{scan_num,proc_type}('Magnetite').area);
            area_Di_grains{scan_num,proc_type}=sum(grains2{scan_num,proc_type}('Diopside').area);
            if max(max(ebsd{scan_num,proc_type}.phase,[],2))>2
            area_Esk_grains{scan_num,proc_type}=sum(grains2{scan_num,proc_type}('Eskolaite').area);
            else
            area_Esk_grains{scan_num,proc_type}=0;
            end
            area_notindexed_grains{scan_num,proc_type}=sum(grains2{scan_num,proc_type}('notIndexed').area);

            sumarea=sum(grains2{scan_num,proc_type}.area);

            afrac_Tmt{scan_num,proc_type}=area_Tmt_grains{scan_num,proc_type}/sumarea;
            afrac_Di{scan_num,proc_type}=area_Di_grains{scan_num,proc_type}/sumarea;
            afrac_Esk{scan_num,proc_type}=area_Esk_grains{scan_num,proc_type}/sumarea;
            afrac_notindexed{scan_num,proc_type}=area_notindexed_grains{scan_num,proc_type}/sumarea;

            % % % Grain sizes (from ellipsoids) % % % %

            % Tmt size information

            [~,majoraxM{scan_num,proc_type},minoraxM{scan_num,proc_type}] = grains2{scan_num,proc_type}('Magnetite').fitEllipse;
            
            % MTEX provides the semi-major and semi-minor axes of the ellipse
            % -> these values are multiplied by 2 to obtain total length
            % and width of ellipse, more commonly used as "grain size" in
            % petrology etc.
            majoraxM{scan_num,proc_type}=majoraxM{scan_num,proc_type}.*2;
            minoraxM{scan_num,proc_type}=minoraxM{scan_num,proc_type}.*2;

            Tmt_areas=grains2{scan_num,proc_type}('Magnetite').area;

            [~,sareaID]=sort(Tmt_areas);
            majors= majoraxM{scan_num,proc_type};
            minors= minoraxM{scan_num,proc_type};
            smaj=majors(sareaID);
            smin=minors(sareaID);

            majax10_M{scan_num,proc_type}=smaj(end-9:end);
            minax10_M{scan_num,proc_type}=smin(end-9:end);
            [std_majax10_M{scan_num,proc_type} , mean_majax10_M{scan_num,proc_type}] = std(smaj(end-9:end));
            [std_minax10_M{scan_num,proc_type} , mean_minax10_M{scan_num,proc_type}] = std(smin(end-9:end));  

            % Lpont = "length" according to Pontesilli et al. 2019 (sqrt L*W)
            [std_Lpont10_M{scan_num,proc_type} , mean_Lpont10_M{scan_num,proc_type}] = std(sqrt(smaj(end-9:end).*smin(end-9:end)));

            % Gmax in units of cm/s
            Gmax_M{scan_num,proc_type} = ( sqrt( (majax10_M{scan_num,proc_type}/10000) .* (minax10_M{scan_num,proc_type}/10000) ) ) / 3600 ;
            [std_Gmax_M{scan_num,proc_type} , mean_Gmax_M{scan_num,proc_type}] = std(Gmax_M{scan_num,proc_type});

            clear smaj majors sareaID minors smin Tmt_areas

            % Cpx size information

            [~,majoraxD{scan_num,proc_type},minoraxD{scan_num,proc_type}] = grains2{scan_num,proc_type}('Diopside').fitEllipse;
            % MTEX provides the semi-major and semi-minor axes of the ellipse
            % -> these values are multiplied by 2 to obtain total length
            % and width of ellipse, more commonly used as "grain size" in
            % petrology etc.
            majoraxD{scan_num,proc_type}=majoraxD{scan_num,proc_type}.*2;
            minoraxD{scan_num,proc_type}=minoraxD{scan_num,proc_type}.*2;

            Cpx_areas=grains2{scan_num,proc_type}('Diopside').area;

            [~,sareaID]=sort(Cpx_areas);
            majors= majoraxD{scan_num,proc_type};
            minors= minoraxD{scan_num,proc_type};
            smaj=majors(sareaID);
            smin=minors(sareaID);

            majax10_D{scan_num,proc_type}=smaj(end-9:end);
            minax10_D{scan_num,proc_type}=smin(end-9:end);
            [std_majax10_D{scan_num,proc_type} , mean_majax10_D{scan_num,proc_type}] = std(smaj(end-9:end));
            [std_minax10_D{scan_num,proc_type} , mean_minax10_D{scan_num,proc_type}] = std(smin(end-9:end));  

            % Lpont = "length" according to Pontesilli et al. 2019 (sqrt L*W)
            [std_Lpont10_D{scan_num,proc_type} , mean_Lpont10_D{scan_num,proc_type}] = std(sqrt(smaj(end-9:end).*smin(end-9:end)));

            % Gmax in units of cm/s
            Gmax_D{scan_num,proc_type} = ( sqrt( (majax10_D{scan_num,proc_type}/10000) .* (minax10_D{scan_num,proc_type}/10000) ) ) / 3600 ;
            [std_Gmax_D{scan_num,proc_type} , mean_Gmax_D{scan_num,proc_type}] = std(Gmax_D{scan_num,proc_type});

            clear smaj majors sareaID minors smin Cpx_areas

            % Boundary lengths and SvP calculation

            GBlength_Di{scan_num,proc_type}=sum(grains2{scan_num,proc_type}('Diopside').boundary.segLength);
            GBlength_Mt{scan_num,proc_type}=sum(grains2{scan_num,proc_type}('Magnetite').boundary.segLength);

            % Cpx SvP in units of mm^-1
            SVPcpx{scan_num,proc_type}= ((4/pi())*(GBlength_Di{scan_num,proc_type}/1000) / ...
            (true_area{scan_num,proc_type}/1e6)) / afrac_Di{scan_num,proc_type} ;

            % Characteristic Tmt size in µm 
            d_Tmt{scan_num,proc_type}= 1000*sqrt( afrac_Tmt{scan_num,proc_type} / N_density_Mt_true{scan_num,proc_type} );
            
            % Gbatch (Tmt) in cm/s
            Gbatch_Tmt{scan_num,proc_type}=d_Tmt{scan_num,proc_type}/(3600*1e5) ;

            clear sumarea x2 y2 xnorm ynorm dunorm dvnorm unorm vnorm unew vnew topleft_id topright_id botright_id botleft_id eyedees ext          
        
    end

end

% Create a scan_num x proc_num sized data structure containing the
% following data, under fieldnames as defined here
datastruct=struct('mapname',mapname,'stepsize',stepsize,'stepsize_regrid',stepsize_regrid,'xspan',xspan,'yspan',yspan,'true_area',true_area,...
    'ebsd',ebsd,'grains2',grains2,...
    'majoraxM',majoraxM,'majoraxD',majoraxD,'minoraxM',minoraxM,'minoraxD',minoraxD,...
    'majax10_M',majax10_M,'minax10_M',minax10_M,'std_majax10_M',std_majax10_M,'mean_majax10_M',mean_majax10_M,...
    'std_minax10_M',std_minax10_M,'mean_minax10_M',mean_minax10_M,...
    'std_Lpont10_M',std_Lpont10_M,'mean_Lpont10_M', mean_Lpont10_M,...
    'Gmax_M',Gmax_M,'std_Gmax_M',std_Gmax_M,'mean_Gmax_M',mean_Gmax_M,...
    'majax10_D',majax10_D,'minax10_D',minax10_D,'std_majax10_D',std_majax10_D,'mean_majax10_D',mean_majax10_D,...
    'std_minax10_D',std_minax10_D,'mean_minax10_D',mean_minax10_D,...
    'std_Lpont10_D',std_Lpont10_D,'mean_Lpont10_D', mean_Lpont10_D,...
    'Gmax_D',Gmax_D,'std_Gmax_D',std_Gmax_D,'mean_Gmax_D',mean_Gmax_D,...
    'NAtmt',N_density_Mt_true,'N_Mt_grains',N_Mt_grains,...
    'area_Tmt_grains',area_Tmt_grains,'area_Di_grains',area_Di_grains,'area_Esk_grains',area_Esk_grains,'area_notindexed_grains',area_notindexed_grains,...
    'afrac_Tmt',afrac_Tmt,'afrac_Di',afrac_Di,'afrac_Esk',afrac_Esk,'afrac_notindexed',afrac_notindexed,...
    'GBlength_Mt',GBlength_Mt,'GBlength_Di',GBlength_Di,...
    'SVPcpx',SVPcpx,'d_Tmt',d_Tmt,'Gbatch_Tmt',Gbatch_Tmt );