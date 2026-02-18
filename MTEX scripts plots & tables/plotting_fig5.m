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
input1=load('ETK5a5b_allscans_distort_regrid.mat','datastruct');

datastruct=input1.datastruct;

clear input1 

fprintf('MATLAB variables loaded!\n\n') 

%% Plotting format settings

lw=5; % linewidth
fs=56; % fontsize x and y label
ax_fs = 46; % fontsize axis numbering
ax_lw = 5; % axis linewidth
msz=28; % markersize

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

% sel variable indicates which results should be plotted (i.e. only plot
% results for step sizes that are equal to or greater than that of the
% original map)
sel={1:6;2:6;2:6;6};

%% Fig 5 a

ylab='\phi_C_p_x (%)';

figure,

    % Loop over each scan and processing type
    for scan_num=1:4
        for proc_type=1:2      
        
            step_regrid(proc_type,:)=[datastruct(scan_num,proc_type).data.stepsize_regrid];

            di_g(proc_type,:)=[datastruct(scan_num,proc_type).data.area_Di_grains];
            tmt_g(proc_type,:)=[datastruct(scan_num,proc_type).data.area_Tmt_grains];
            esk_g(proc_type,:)=[datastruct(scan_num,proc_type).data.area_Esk_grains];
            not_g(proc_type,:)=[datastruct(scan_num,proc_type).data.area_notindexed_grains];
        
            sum_g(proc_type,:)=di_g(proc_type,:)+tmt_g(proc_type,:)+esk_g(proc_type,:)+not_g(proc_type,:);
            
            % %
            afrac_di_g(proc_type,:)=di_g(proc_type,:)./sum_g(proc_type,:);
            afrac_tmt_g(proc_type,:)=tmt_g(proc_type,:)./sum_g(proc_type,:);
            afrac_esk_g(proc_type,:)=esk_g(proc_type,:)./sum_g(proc_type,:);
            afrac_not_g(proc_type,:)=not_g(proc_type,:)./sum_g(proc_type,:);
        
        end
        
        % two lines are plotted, one for standardised, one for
        % non-standardised scan

        plot(step_regrid(2,sel{scan_num}),100*(afrac_di_g(2,sel{scan_num})),'marker','o','Color',color(scan_num,:),'MarkerFaceColor','w','LineStyle',':','LineWidth',lw,'MarkerSize',msz)    
        hold on
        plot(step_regrid(1,sel{scan_num}),100*(afrac_di_g(1,sel{scan_num})),'marker','x','Color',color(scan_num,:),'LineWidth',lw,'MarkerSize',msz)
        
        xlim([0.1 1.1])
        xticks([0.15 0.2 0.4 0.6 0.8 1])
        xlabel('step size (µm)')
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

exportgraphics(f,'Cpx_crystallinity_regrid.png','Resolution',300)

%% Fig 2 b

ylab='\phi_T_m_t (%)';

figure,

    % Loop over each scan and processing type
    for scan_num=1:4
        for proc_type=1:2    

            step_regrid(proc_type,:)=[datastruct(scan_num,proc_type).data.stepsize_regrid];
        
            di_g(proc_type,:)=[datastruct(scan_num,proc_type).data.area_Di_grains];
            tmt_g(proc_type,:)=[datastruct(scan_num,proc_type).data.area_Tmt_grains];
            esk_g(proc_type,:)=[datastruct(scan_num,proc_type).data.area_Esk_grains];
            not_g(proc_type,:)=[datastruct(scan_num,proc_type).data.area_notindexed_grains];
        
            sum_g(proc_type,:)=di_g(proc_type,:)+tmt_g(proc_type,:)+esk_g(proc_type,:)+not_g(proc_type,:);
            
            % %
            afrac_di_g(proc_type,:)=di_g(proc_type,:)./sum_g(proc_type,:);
            afrac_tmt_g(proc_type,:)=tmt_g(proc_type,:)./sum_g(proc_type,:);
            afrac_esk_g(proc_type,:)=esk_g(proc_type,:)./sum_g(proc_type,:);
            afrac_not_g(proc_type,:)=not_g(proc_type,:)./sum_g(proc_type,:);
        
        end
        
        % two lines are plotted, one for standardised, one for
        % non-standardised scan

        plot(step_regrid(2,sel{scan_num}),100*(afrac_tmt_g(2,sel{scan_num})),'marker','o','Color',color(scan_num,:),'MarkerFaceColor','w','LineStyle',':','LineWidth',lw,'MarkerSize',msz)    
        hold on
        plot(step_regrid(1,sel{scan_num}),100*(afrac_tmt_g(1,sel{scan_num})),'marker','x','Color',color(scan_num,:),'LineWidth',lw,'MarkerSize',msz)
        
        xlim([0.1 1.1])
        xticks([0.15 0.2 0.4 0.6 0.8 1])
        xlabel('step size (µm)')
        ylim([0.5 3])
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

exportgraphics(f,'Tmt_crystallinity_regrid.png','Resolution',300)


%% Fig 2 d

ylab='L_m_a_x Cpx (\mum)';

figure,

    % Loop over each scan and processing type
    for scan_num=1:4
        for proc_type=1:2      

            step_regrid(proc_type,:)=[datastruct(scan_num,proc_type).data.stepsize_regrid];
        
            for ii=1:length(step_regrid(proc_type,:))
                if length(datastruct(scan_num,proc_type).data(ii).regrid_grains) > 1 % needed as for "dummy" data length=1
                majoraxD=datastruct(scan_num,proc_type).data(ii).majoraxD;
                [~,sort_areaID]=sort(datastruct(scan_num,proc_type).data(ii).regrid_grains('Diopside').area);
                majoraxD=majoraxD(sort_areaID);
                lmax10_D(ii)=mean([majoraxD(end-9:end)]);
                else
                lmax10_D(ii)=1 ; 
                end
            end
            means(proc_type,:)=lmax10_D;
            
        end
        
        % two lines are plotted, one for standardised, one for
        % non-standardised scan

        plot(step_regrid(2,sel{scan_num}),means(2,sel{scan_num}),'marker','o','Color',color(scan_num,:),'MarkerFaceColor','w','LineStyle',':','LineWidth',lw,'MarkerSize',msz)    
        hold on
        plot(step_regrid(1,sel{scan_num}),means(1,sel{scan_num}),'marker','x','Color',color(scan_num,:),'LineWidth',lw,'MarkerSize',msz)
        
        xlim([0.1 1.1])
        xticks([0.15 0.2 0.4 0.6 0.8 1])
        xlabel('step size (µm)')
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

exportgraphics(f,'Cpx_Lmax_regrid.png','Resolution',300)

%% Fig 2 d

ylab='L_m_a_x Tmt (\mum)';

figure,

    % Loop over each scan and processing type
    for scan_num=1:4
        for proc_type=1:2  

            step_regrid(proc_type,:)=[datastruct(scan_num,proc_type).data.stepsize_regrid];
        
            for ii=1:length(step_regrid(proc_type,:))
                if length(datastruct(scan_num,proc_type).data(ii).regrid_grains) > 1 % needed as for "dummy" data length=1
                majoraxM=datastruct(scan_num,proc_type).data(ii).majoraxM;
                [~,sort_areaID]=sort(datastruct(scan_num,proc_type).data(ii).regrid_grains('Magnetite').area);
                majoraxM=majoraxM(sort_areaID);
                lmax10_M(ii)=mean([majoraxM(end-9:end)]);
                else
                lmax10_M(ii)=1 ; 
                end
            end    
            means(proc_type,:)=lmax10_M;
            
        end
        
        % two lines are plotted, one for standardised, one for
        % non-standardised scan

        plot(step_regrid(2,sel{scan_num}),means(2,sel{scan_num}),'marker','o','Color',color(scan_num,:),'MarkerFaceColor','w','LineStyle',':','LineWidth',lw,'MarkerSize',msz)    
        hold on
        plot(step_regrid(2,sel{scan_num}),means(1,sel{scan_num}),'marker','x','Color',color(scan_num,:),'LineWidth',lw,'MarkerSize',msz)
        
        xlim([0.1 1.1])
        xticks([0.15 0.2 0.4 0.6 0.8 1])
        xlabel('step size (µm)')
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

exportgraphics(f,'Tmt_Lmax_regrid.png','Resolution',300)

%% Fig 2 b

ylab='S_v^P Cpx (mm^-^1)';

figure,

    % Loop over each scan and processing type
    for scan_num=1:4
        for proc_type=1:2

            step_regrid(proc_type,:)=[datastruct(scan_num,proc_type).data.stepsize_regrid];
        
            gbldi(proc_type,:)=[datastruct(scan_num,proc_type).data.GBlength_Di];
        
            di_g(proc_type,:)=[datastruct(scan_num,proc_type).data.area_Di_grains];
            tmt_g(proc_type,:)=[datastruct(scan_num,proc_type).data.area_Tmt_grains];
            esk_g(proc_type,:)=[datastruct(scan_num,proc_type).data.area_Esk_grains];
            not_g(proc_type,:)=[datastruct(scan_num,proc_type).data.area_notindexed_grains];
            
            sum_g(proc_type,:)=di_g(proc_type,:)+tmt_g(proc_type,:)+esk_g(proc_type,:)+not_g(proc_type,:);
            
            % %
            afrac_di_g(proc_type,:)=di_g(proc_type,:)./sum_g(proc_type,:);
            afrac_tmt_g(proc_type,:)=tmt_g(proc_type,:)./sum_g(proc_type,:);
            afrac_esk_g(proc_type,:)=esk_g(proc_type,:)./sum_g(proc_type,:);
            afrac_not_g=not_g(proc_type,:)./sum_g(proc_type,:);
        
            svp(proc_type,:)=((4/pi())*(gbldi(proc_type,:)/1000)./...
            (datastruct(scan_num).true_area/1e6))./afrac_di_g(proc_type,:); 
        
        end
        
        % two lines are plotted, one for standardised, one for
        % non-standardised scan

        plot(step_regrid(2,sel{scan_num}),svp(2,sel{scan_num}),'marker','o','Color',color(scan_num,:),'MarkerFaceColor','w','LineStyle',':','LineWidth',lw,'MarkerSize',msz)    
        hold on
        plot(step_regrid(1,sel{scan_num}),svp(1,sel{scan_num}),'marker','x','Color',color(scan_num,:),'LineWidth',lw,'MarkerSize',msz)
        
        xlim([0.1 1.1])
        xticks([0.15 0.2 0.4 0.6 0.8 1])
        xlabel('step size (µm)')
        ylim([300 1300])
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

exportgraphics(f,'Cpx_SvP_regrid.png','Resolution',300)

%% Fig 2 b

ylab='N_A Tmt (mm^-^2)';

figure,

    % Loop over each scan and processing type
    for scan_num=1:4
        for proc_type=1:2

            step_regrid(proc_type,:)=[datastruct(scan_num,proc_type).data.stepsize_regrid];
        
            ndens_tmt(proc_type,:)=[datastruct(scan_num,proc_type).data.N_density_Mt_true];
        
        end
        
        % two lines are plotted, one for standardised, one for
        % non-standardised scan

        plot(step_regrid(2,sel{scan_num}),ndens_tmt(2,sel{scan_num}),'marker','o','Color',color(scan_num,:),'MarkerFaceColor','w','LineStyle',':','LineWidth',lw,'MarkerSize',msz)    
        hold on
        plot(step_regrid(1,sel{scan_num}),ndens_tmt(1,sel{scan_num}),'marker','x','Color',color(scan_num,:),'LineWidth',lw,'MarkerSize',msz)
        
        xlim([0.1 1.1])
        xticks([0.15 0.2 0.4 0.6 0.8 1])
        xlabel('step size (µm)')
        ylim([1000 5000])
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

exportgraphics(f,'Tmt_NA_regrid.png','Resolution',300)