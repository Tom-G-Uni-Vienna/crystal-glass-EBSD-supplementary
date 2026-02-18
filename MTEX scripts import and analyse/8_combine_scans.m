%% Specify Plotting Convention and Load data

% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

% Enter the location of the supplementary data folder
supplementary_location =...
    'Z:\Crystal clustering project RESEARCH\Papers\EBSD vs BSE paper';

% Output path (automatically derived as long as supp. mat. structure
% remains unchanged)
savelocation = [ supplementary_location...
    '\FINAL Supplementary\FINAL scripts & data\MATLAB data'];

% Path to load MATLAB variables (automatically derived as long as supp. mat. structure
% remains unchanged)
loadlocation = savelocation;

fprintf('\nLoading MATLAB variables...\n') 

cd(loadlocation)
% MATLAB prefers to load data to structures and extract variables from them
% instead of loading the variables directly 
% (https://de.mathworks.com/matlabcentral/answers/28676-why-use-x-load-myfile-mat)
input1=load('ETK5a5b_allscans_distort_regrid_final.mat','datastruct');
input2=load('ETK5a5b_allscans_distort_regrid_final_075.mat','datastruct075');
input3=load('ETK5a5b_allscans_distort_regrid_final_125.mat','datastruct125');

datastruct=input1.datastruct;
datastruct075=input2.datastruct075;
datastruct125=input3.datastruct125;

clear input1 input2 input3

fprintf('MATLAB variables loaded!\n\n') 


%% Combining data from multiple maps on the same sample

fprintf('\ncombining map results...\n') 

% Create an empty structure to output results into - 2 samples x 2
% processing types (rows = samples, columns = processing types)
resultstruct(2,2)=struct;

% Helper variable for the combination of data - needed as 
% sample 1 = maps 1+2, sample 2 = maps 3+4 ->
% so the two maps for a given sample in the input variable datastruct are
% at row positions q(sample_num) and q(sample_num)+1.
q=[1,3]; 

% % % Part 1: combining the variables listed in res_list below % % % % 
res_list={'true_area','N_Mt_grains',...
    'area_Tmt_grains','area_Di_grains','area_Esk_grains','area_notindexed_grains',...
    'GBlength_Mt','GBlength_Di'};

for sample_num=1:2
    for proc_type=1:2
        % Cycle over all variables in res_list
        for phase_num=1:length(res_list)
            % For each position in resultstruct, and for each variable,
            % create THREE entries - number (1) is the combined result for CI
            % threshold = 0.1, number (2) is the combined result for CI threshold
            % 0.075, and (3) is the combined result for CI threshold 0.125.
            % e.g. to later access the combined true_area result for sample 1, CI
            % standardised, CI threshold 0.125 we would write:
            % resultstruct(1,2).true_area(3)
            resultstruct(sample_num,proc_type).(res_list{phase_num})(1) = sum([ datastruct(q(sample_num),proc_type).(res_list{phase_num}); datastruct(q(sample_num)+1,proc_type).(res_list{phase_num}) ]);
            resultstruct(sample_num,proc_type).(res_list{phase_num})(2) = sum([ datastruct075(q(sample_num),proc_type).(res_list{phase_num}); datastruct075(q(sample_num)+1,proc_type).(res_list{phase_num}) ]);
            resultstruct(sample_num,proc_type).(res_list{phase_num})(3) = sum([ datastruct125(q(sample_num),proc_type).(res_list{phase_num}); datastruct125(q(sample_num)+1,proc_type).(res_list{phase_num}) ]);
        end
    end
end

% % % Part 2: calculating overall NA_tmt, SvP_cpx, dTmt, Gbatch % % % % 

% For these we need individual calculations for each sample and processing type
% using the sums generated above.
% N.B. NA_tmt etc. uses true_area/1e6 to convert to mm^2

phaselist={'Di','Tmt','Esk','notindexed'};

for sample_num=1:2
    for proc_type=1:2 

        % Create temp variable containing the data to be processed in this
        % step - not strictly necessary but makes code easier to read
        % below...
        res=resultstruct(sample_num,proc_type);

        % Cycle over phase list above, for each phase calculate area
        % fraction: i.e. grain area of that phase divided by the sum of 
        % all phase areas, append this info to resultstruct
        for i=1:length(phaselist)
            resultstruct(sample_num,proc_type).(['afrac_' phaselist{i}]) = ...
                res.(['area_' phaselist{i} '_grains'])./ ...
                    (res.area_Di_grains + ...
                    res.area_Tmt_grains + ...
                    res.area_Esk_grains + ...
                    res.area_notindexed_grains) ;
        end 

        % Add NA_tmt to resultstruct (Total Ntmt / total map area)
        resultstruct(sample_num,proc_type).NAtmt=res.N_Mt_grains./(res.true_area/1e6) ;

        % Add SvPcpx to resultstruct (function of total cpx boundary length
        % and total map area)
        resultstruct(sample_num,proc_type).SVPcpx = ...
            ((4/pi()).*(res.GBlength_Di/1000)./ ...
            (res.true_area/1e6))./resultstruct(sample_num,proc_type).afrac_Di ;

        % Add d_tmt to resultstruct (function of overall tmt area fraction 
        % and NA_tmt)
        resultstruct(sample_num,proc_type).d_Tmt =...
            1000*sqrt( resultstruct(sample_num,proc_type).afrac_Tmt ./...
            resultstruct(sample_num,proc_type).NAtmt );
        
        % Add Gbatch_tmt to resultstruct (function of d_tmt and experimental
        % duration, factor 1e5 to convert to cms^-1).
        resultstruct(sample_num,proc_type).Gbatch_Tmt = resultstruct(sample_num,proc_type).d_Tmt/(3600*1e5);

        clear res

    end
end

% % % Part 3: calculating overall Lmax10, sqrt(L*W)max10, Gmax10, and 
% standard deviations of these parameters % % % % % % % % % % % 

p_list={'M','D'};
p_list2={'Magnetite','Diopside'};

% Looping over each sample and each processing type...
for sample_num=1:2
    for proc_type=1:2

        % Create two temporary structure arrays d1 and d2, each with
        % three entries. d1 contains all data for the first map on a given
        % sample, d2 all data for the second map, the three entries 
        % correspond to data for CI thresholds 0.1, 0.075, and 0.125, 
        % respectively. Not strictly necessary, but enhances the readability
        % of the following code. 
        d1(1)=datastruct(q(sample_num),proc_type);
        d1(2)=datastruct075(q(sample_num),proc_type);
        d1(3)=datastruct125(q(sample_num),proc_type);
    
        d2(1)=datastruct(q(sample_num)+1,proc_type);
        d2(2)=datastruct075(q(sample_num)+1,proc_type);
        d2(3)=datastruct125(q(sample_num)+1,proc_type);       

        % Repeat the code below first for Tmt ('M') and then for Cpx ('D')
        for phase_num=1:length(p_list)            

            % Cycle through the three sets of results (for the three CI
            % thresholds)
            for ci_val=1:3    
 
                % Combine all major and minor axes for both scans on a
                % sample
                majors=[ d1(ci_val).(['majorax' p_list{phase_num}]); d2(ci_val).(['majorax' p_list{phase_num}]) ];
                minors=[ d1(ci_val).(['minorax' p_list{phase_num}]); d2(ci_val).(['minorax' p_list{phase_num}]) ];

                % Combine grain areas for both scans on a sample
                areas=[ d1(ci_val).grains2(p_list2{phase_num}).area ; d2(ci_val).grains2(p_list2{phase_num}).area ];

                % Determine ranking of grains by area
                [~,sareasID]=sort(areas);

                % Sort minor and major axes by the area of the grain they
                % belong to
                smaj=majors(sareasID);
                smin=minors(sareasID);

                % Calculate the sqrt(L*W) for all grains (sorted by area)
                Lpont_byArea=sqrt(smaj.*smin);
               
                % Add the combined major, minor, and area data to the
                % result structure
                resultstruct(sample_num,proc_type).(['majors_' p_list{phase_num}]){ci_val}=majors;
                resultstruct(sample_num,proc_type).(['minors_' p_list{phase_num}]){ci_val}=minors;
                resultstruct(sample_num,proc_type).(['areas_' p_list{phase_num}]){ci_val}=areas;    
        
                % Add the major and minor axes of the largest 10 grains by
                % area (considering both scans on the sample!) to the results structure
                resultstruct(sample_num,proc_type).(['majax10_' p_list{phase_num}]){ci_val}=smaj(end-9:end);
                resultstruct(sample_num,proc_type).(['minax10_' p_list{phase_num}]){ci_val}=smin(end-9:end);
      
                % Calculate standard deviation and mean of the major axes 
                % of the 10 largest crystals (using command std, which outputs stdev and mean) 
                [resultstruct(sample_num,proc_type).(['std_majax10_' p_list{phase_num}])(ci_val),...
                    resultstruct(sample_num,proc_type).(['mean_majax10_' p_list{phase_num}])(ci_val) ] =...
                    std(smaj(end-9:end));
        
                % Calculate standard deviation and mean of the minor axes 
                % of the 10 largest crystals (using command std, which outputs stdev and mean) 
                [resultstruct(sample_num,proc_type).(['std_minax10_' p_list{phase_num}])(ci_val),...
                    resultstruct(sample_num,proc_type).(['mean_minax10_' p_list{phase_num}])(ci_val) ] =...
                    std(smin(end-9:end));  

                % Calculate standard deviation and mean of the sqrt(L*W) 
                % of the 10 largest crystals (using command std, which outputs stdev and mean) 
                [resultstruct(sample_num,proc_type).(['std_L10pont_Area_' p_list{phase_num}])(ci_val),...
                    resultstruct(sample_num,proc_type).(['mean_L10pont_Area_' p_list{phase_num}])(ci_val) ] =...
                    std(Lpont_byArea(end-9:end)); 
        
                % Calculate the growth rates of the 10 largest crystals 
                % (by area), factor 1e5 converts to cm/s 
                % (3600 s = experimental duration of 30 mins)
                resultstruct(sample_num,proc_type).(['Gmax_' p_list{phase_num}]){ci_val} = ...
                    ( sqrt( (resultstruct(sample_num,proc_type).(['majax10_' p_list{phase_num}]){ci_val}/10000) ...
                    .* (resultstruct(sample_num,proc_type).(['minax10_' p_list{phase_num}]){ci_val}/10000)  ) ) / 3600 ;

                % Calculate standard deviation and mean of the growth rates
                % of the 10 largest crystals (using command std, which outputs stdev and mean) 
                [resultstruct(sample_num,proc_type).(['stdGmax_' p_list{phase_num}])(ci_val),...
                    resultstruct(sample_num,proc_type).(['meanGmax_' p_list{phase_num}])(ci_val)] = ...
                    std(resultstruct(sample_num,proc_type).(['Gmax_' p_list{phase_num}]){ci_val});
        
                clear majors minors smin smaj areas Lpont_byArea sareasID

            end

        end
    end
end

clear d1 d2  datastruct datastruct125 datastruct075

fprintf('\nsaving results...\n') 

cd(savelocation)

save('ETK5a5b_allscans_distort_combined_result.mat','resultstruct')

fprintf('done!\n\n') 

