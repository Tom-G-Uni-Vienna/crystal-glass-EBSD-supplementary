%% Load data & setup formatting

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

fprintf('loading MATLAB variables...\n\n') 

cd(loadlocation)

% MATLAB prefers to load data to structures and extract variables from them
% instead of loading the variables directly 
% (https://de.mathworks.com/matlabcentral/answers/28676-why-use-x-load-myfile-mat)
input1=load('ETK5a5b_allscans_distort_grainsmooth.mat','datastruct');

datastruct=input1.datastruct;

clear input1 

fprintf('MATLAB variables loaded!\n\n') 

smoothnum=[datastruct(1).smooth_result.smooth_num];

%% Plotting format settings

lw=4.5; % linewidth
fs=43; % fontsize x and y label
ax_fs = 35; % fontsize axis numbering
ax_lw = 4.5; % axis linewidth
msz=20; % markersize

% Get screen sizes and define screen fractions for plot positioning
screendat=get(0,'ScreenSize');
w=screendat(3); % width
h=screendat(4); % height
leftfrac=0.01;
botfrac=0.02;
horzfrac=0.6;
vertfrac=0.86;

vert_div_horiz=1;

% Generate color list of desired order
color=colororder('gem');
c2=color(2,:);
c4=color(4,:);
color(4,:)=c2;
color(2,:)=c4;
color=flip(color);
color=color([4:7],:);



%% Fig 2 f

ylab='S_v^P Cpx (mm^-^1)';

figure,

    % Loop over each scan and processing type
    for scan_num=1:4
        for proc_type=1:2      
        
         gbldi(proc_type,:)=[datastruct(scan_num,proc_type).smooth_result.GBlength_Di];

        di_g(proc_type,:)=[datastruct(scan_num,proc_type).smooth_result.area_Di_grains];
        tmt_g(proc_type,:)=[datastruct(scan_num,proc_type).smooth_result.area_Tmt_grains];
        esk_g(proc_type,:)=[datastruct(scan_num,proc_type).smooth_result.area_Esk_grains];
        not_g(proc_type,:)=[datastruct(scan_num,proc_type).smooth_result.area_notindexed_grains];

        sum_g(proc_type,:)=di_g(proc_type,:)+tmt_g(proc_type,:)+esk_g(proc_type,:)+not_g(proc_type,:);

        % %
        afrac_di_g(proc_type,:)=di_g(proc_type,:)./sum_g(proc_type,:);
        afrac_tmt_g(proc_type,:)=tmt_g(proc_type,:)./sum_g(proc_type,:);
        afrac_esk_g(proc_type,:)=esk_g(proc_type,:)./sum_g(proc_type,:);
        afrac_not_g=not_g(proc_type,:)./sum_g(proc_type,:);

        svp(proc_type,:)=((4/pi())*(gbldi(proc_type,:)/1000)./...
            (datastruct(scan_num).true_area/1e6))./afrac_di_g(proc_type,:);
            
        end
        
        % jj gives column number of the Fit or CI threshold matrix to be
        % plotted, and selects corresponding values to be plotted on y-axis
        % two lines are plotted, one for standardised, one for
        % non-standardised scan

        plot(smoothnum,svp(2,:),'marker','o','Color',color(scan_num,:),'MarkerFaceColor','w','LineStyle',':','LineWidth',lw,'MarkerSize',msz)    
        hold on
        plot(smoothnum,svp(1,:),'marker','x','Color',color(scan_num,:),'LineWidth',lw,'MarkerSize',msz)
        
        xlim([0 21])
        xticks(smoothnum)
        xlabel('N smoothing iterations')
        ylim([300 1700])
        ylabel(ylab)
        
        ax = gca;
        ax.XAxis.FontSize = ax_fs;
        ax.YAxis.FontSize = ax_fs;
        ax.XLabel.FontSize = fs-3;
        % ax.XLabel.FontWeight = 'bold';
        % ax.XLabel.Position(2) = ax.XLabel.Position(2)+0.25
        ax.XLabel.FontAngle = 'italic';
        ax.XTickLabelRotation = 90;
        ax.YLabel.FontSize = fs;
        % ax.YLabel.FontWeight = 'bold';
        % ax.YLabel.FontAngle = 'italic';
        ax.LineWidth = ax_lw;
     end

cd(savelocation)

f = gcf;
f.Position = [leftfrac*w botfrac*h horzfrac*w horzfrac*vert_div_horiz*w];

exportgraphics(f,'Cpx_SvP_smooths.png','Resolution',300)