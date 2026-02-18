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

lw=5; % linewidth
fs=53; % fontsize x and y label
ax_fs = 38; % fontsize axis numbering
ax_lw = 5; % axis linewidth
msz=22; % markersize

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

%% Fig 2 a

ylab='\phi_C_p_x (%)';

figure,

    % Loop over each scan and processing type
    for scan_num=1:4
        for proc_type=1:2      
        
            % For each scan & proc type, obtain a matrix of results the
            % same size as the matrix of input smoothnumthresh,
            % store these all together in one indexable matrix
            di_g(proc_type,:)=[datastruct(scan_num,proc_type).smooth_result.area_Di_grains];
            tmt_g(proc_type,:)=[datastruct(scan_num,proc_type).smooth_result.area_Tmt_grains];
            esk_g(proc_type,:)=[datastruct(scan_num,proc_type).smooth_result.area_Esk_grains];
            not_g(proc_type,:)=[datastruct(scan_num,proc_type).smooth_result.area_notindexed_grains];
        
            sum_g(proc_type,:)=di_g(proc_type,:)+tmt_g(proc_type,:)+esk_g(proc_type,:)+not_g(proc_type,:);
            
            % %
            afrac_di_g(proc_type,:)=di_g(proc_type,:)./sum_g(proc_type,:);
            afrac_tmt_g(proc_type,:)=tmt_g(proc_type,:)./sum_g(proc_type,:);
            afrac_esk_g(proc_type,:)=esk_g(proc_type,:)./sum_g(proc_type,:);
            afrac_not_g(proc_type,:)=not_g(proc_type,:)./sum_g(proc_type,:);
        
        end
        
        % jj gives column number of the Fit or CI threshold matrix to be
        % plotted, and selects corresponding values to be plotted on y-axis
        % two lines are plotted, one for standardised, one for
        % non-standardised scan

        plot(smoothnum,100*(afrac_di_g(2,:)),'marker','o','Color',color(scan_num,:),'MarkerFaceColor','w','LineStyle',':','LineWidth',lw,'MarkerSize',msz)    
        hold on
        plot(smoothnum,100*(afrac_di_g(1,:)),'marker','x','Color',color(scan_num,:),'LineWidth',lw,'MarkerSize',msz)
        
        xlim([0 21])
        xticks(smoothnum)
        xlabel('N smoothing iterations')
        ylim([30 70])
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

exportgraphics(f,'Cpx_crystallinity_smooths.png','Resolution',300)

%% Fig 2 b

ylab='\phi_T_m_t (%)';

figure,

    % Loop over each scan and processing type
    for scan_num=1:4
        for proc_type=1:2      
        
            % For each scan & proc type, obtain a matrix of results the
            % same size as the matrix of input smoothnumthresh,
            % store these all together in one indexable matrix
            di_g(proc_type,:)=[datastruct(scan_num,proc_type).smooth_result.area_Di_grains];
            tmt_g(proc_type,:)=[datastruct(scan_num,proc_type).smooth_result.area_Tmt_grains];
            esk_g(proc_type,:)=[datastruct(scan_num,proc_type).smooth_result.area_Esk_grains];
            not_g(proc_type,:)=[datastruct(scan_num,proc_type).smooth_result.area_notindexed_grains];
        
            sum_g(proc_type,:)=di_g(proc_type,:)+tmt_g(proc_type,:)+esk_g(proc_type,:)+not_g(proc_type,:);
            
            % %
            afrac_di_g(proc_type,:)=di_g(proc_type,:)./sum_g(proc_type,:);
            afrac_tmt_g(proc_type,:)=tmt_g(proc_type,:)./sum_g(proc_type,:);
            afrac_esk_g(proc_type,:)=esk_g(proc_type,:)./sum_g(proc_type,:);
            afrac_not_g(proc_type,:)=not_g(proc_type,:)./sum_g(proc_type,:);
        
        end
        
        % jj gives column number of the Fit or CI threshold matrix to be
        % plotted, and selects corresponding values to be plotted on y-axis
        % two lines are plotted, one for standardised, one for
        % non-standardised scan

        plot(smoothnum,100*(afrac_tmt_g(2,:)),'marker','o','Color',color(scan_num,:),'MarkerFaceColor','w','LineStyle',':','LineWidth',lw,'MarkerSize',msz)    
        hold on
        plot(smoothnum,100*(afrac_tmt_g(1,:)),'marker','x','Color',color(scan_num,:),'LineWidth',lw,'MarkerSize',msz)
        
        xlim([0 21])
        xticks(smoothnum)
        xlabel('N smoothing iterations')
        ylim([1 3.5])
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

exportgraphics(f,'Tmt_crystallinity_smooths.png','Resolution',300)

%% Fig 2 c

ylab='\phi_E_s_k (%)';

figure,

    % Loop over each scan and processing type
    for scan_num=1:4
        for proc_type=1:2      
        
            % For each scan & proc type, obtain a matrix of results the
            % same size as the matrix of input smoothnumthresh,
            % store these all together in one indexable matrix
            di_g(proc_type,:)=[datastruct(scan_num,proc_type).smooth_result.area_Di_grains];
            tmt_g(proc_type,:)=[datastruct(scan_num,proc_type).smooth_result.area_Tmt_grains];
            esk_g(proc_type,:)=[datastruct(scan_num,proc_type).smooth_result.area_Esk_grains];
            not_g(proc_type,:)=[datastruct(scan_num,proc_type).smooth_result.area_notindexed_grains];
        
            sum_g(proc_type,:)=di_g(proc_type,:)+tmt_g(proc_type,:)+esk_g(proc_type,:)+not_g(proc_type,:);
            
            % %
            afrac_di_g(proc_type,:)=di_g(proc_type,:)./sum_g(proc_type,:);
            afrac_tmt_g(proc_type,:)=tmt_g(proc_type,:)./sum_g(proc_type,:);
            afrac_esk_g(proc_type,:)=esk_g(proc_type,:)./sum_g(proc_type,:);
            afrac_not_g(proc_type,:)=not_g(proc_type,:)./sum_g(proc_type,:);
        
        end
        
        % jj gives column number of the Fit or CI threshold matrix to be
        % plotted, and selects corresponding values to be plotted on y-axis
        % two lines are plotted, one for standardised, one for
        % non-standardised scan

        plot(smoothnum,100*(afrac_esk_g(2,:)),'marker','o','Color',color(scan_num,:),'MarkerFaceColor','w','LineStyle',':','LineWidth',lw,'MarkerSize',msz)    
        hold on
        plot(smoothnum,100*(afrac_esk_g(1,:)),'marker','x','Color',color(scan_num,:),'LineWidth',lw,'MarkerSize',msz)
        
        xlim([0 21])
        xticks(smoothnum)
        xlabel('N smoothing iterations')
        ylim([0 0.15])
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

exportgraphics(f,'Esk_crystallinity_smooths.png','Resolution',300)

%% Fig 2 d

ylab='L_m_a_x Cpx (\mum)';

figure,

    % Loop over each scan and processing type
    for scan_num=1:4
        for proc_type=1:2      
        
            for ii=1:length(smoothnum)
                majoraxD=datastruct(scan_num,proc_type).smooth_result(ii).majoraxD;
                [~,sort_areaID]=sort(datastruct(scan_num,proc_type).smooth_result(ii).smooth_grains('Diopside').area);
                majoraxD=majoraxD(sort_areaID);
                lmax10_D(ii)=mean([majoraxD(end-9:end)]);
            end
            means(proc_type,:)=lmax10_D;
            
        end
        
        % jj gives column number of the Fit or CI threshold matrix to be
        % plotted, and selects corresponding values to be plotted on y-axis
        % two lines are plotted, one for standardised, one for
        % non-standardised scan

        plot(smoothnum,means(2,:),'marker','o','Color',color(scan_num,:),'MarkerFaceColor','w','LineStyle',':','LineWidth',lw,'MarkerSize',msz)    
        hold on
        plot(smoothnum,means(1,:),'marker','x','Color',color(scan_num,:),'LineWidth',lw,'MarkerSize',msz)
        
        xlim([0 21])
        xticks(smoothnum)
        xlabel('N smoothing iterations')
        ylim([50 130])
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

exportgraphics(f,'Cpx_Lmax_smooths.png','Resolution',300)

%% Fig 2 e

ylab='L_m_a_x Tmt (\mum)';

figure,

    % Loop over each scan and processing type
    for scan_num=1:4
        for proc_type=1:2      
        
            for ii=1:length(smoothnum)
                majoraxM=datastruct(scan_num,proc_type).smooth_result(ii).majoraxM;
                [~,sort_areaID]=sort(datastruct(scan_num,proc_type).smooth_result(ii).smooth_grains('Magnetite').area);
                majoraxM=majoraxM(sort_areaID);
                lmax10_M(ii)=mean([majoraxM(end-9:end)]);
            end
            means(proc_type,:)=lmax10_M;
            
        end
        
        % jj gives column number of the Fit or CI threshold matrix to be
        % plotted, and selects corresponding values to be plotted on y-axis
        % two lines are plotted, one for standardised, one for
        % non-standardised scan

        plot(smoothnum,means(2,:),'marker','o','Color',color(scan_num,:),'MarkerFaceColor','w','LineStyle',':','LineWidth',lw,'MarkerSize',msz)    
        hold on
        plot(smoothnum,means(1,:),'marker','x','Color',color(scan_num,:),'LineWidth',lw,'MarkerSize',msz)
        
        xlim([0 21])
        xticks(smoothnum)
        xlabel('N smoothing iterations')
        ylim([5 13])
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

exportgraphics(f,'Tmt_Lmax_smooths.png','Resolution',300)

%% Fig 2 f

ylab='N_A Tmt (mm^-^2)';

figure,

    % Loop over each scan and processing type
    for scan_num=1:4
        for proc_type=1:2      
        
            ndens_tmt(proc_type,:)=[datastruct(scan_num,proc_type).smooth_result.N_density_Mt_true];
            
        end
        
        % jj gives column number of the Fit or CI threshold matrix to be
        % plotted, and selects corresponding values to be plotted on y-axis
        % two lines are plotted, one for standardised, one for
        % non-standardised scan

        plot(smoothnum,ndens_tmt(2,:),'marker','o','Color',color(scan_num,:),'MarkerFaceColor','w','LineStyle',':','LineWidth',lw,'MarkerSize',msz)    
        hold on
        plot(smoothnum,ndens_tmt(1,:),'marker','x','Color',color(scan_num,:),'LineWidth',lw,'MarkerSize',msz)
        
        xlim([0 21])
        xticks(smoothnum)
        xlabel('N smoothing iterations')
        ylim([500 5500])
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

exportgraphics(f,'Tmt_NA_smooths.png','Resolution',300)