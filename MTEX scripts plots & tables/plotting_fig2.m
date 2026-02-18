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
input1=load('ETK5a5b_allscans_distort_cithresh.mat','datastruct');

datastruct=input1.datastruct;

clear input1 

fprintf('MATLAB variables loaded!\n\n') 

syze=size(datastruct(1).cifit_result);

cithresh=datastruct(1).cithresh;
fitthresh=datastruct(1).fitthresh;

%% Plotting format settings

lw=2.5; % linewidth
fs=28; % fontsize x and y label
ax_fs = 23; % fontsize axis numbering
ax_lw = 2.5; % axis linewidth
msz=16; % markersize

% Get screen sizes and define screen fractions for plot positioning
screendat=get(0,'ScreenSize');
w=screendat(3); % width
h=screendat(4); % height
leftfrac=0.01;
botfrac=0.02;
horzfrac=0.98;
vertfrac=0.86;

vert_div_horiz=0.37

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

% titlelist={'CI thresh only','Fit thresh only','CI thresh (Fit=1.4)','Fit Thresh (CI=0.1)'};

figure,
% For 4 subplots...
t=tiledlayout(1,4);

for jj=1:4

nexttile    
% subplot(1,4,jj)

    % Loop over each scan and processing type
    for scan_num=1:4
        for proc_type=1:2      
        
            % For each scan & proc type, obtain a 5x4 matrix of results the
            % same size as the matrix of input CI and Fit combinations,
            % store these all together in one indexable 4D matrix
            di_g(:,:,proc_type,scan_num)=reshape([datastruct(scan_num,proc_type).cifit_result.area_Di_grains],syze);
            tmt_g(:,:,proc_type,scan_num)=reshape([datastruct(scan_num,proc_type).cifit_result.area_Tmt_grains],syze);
            esk_g(:,:,proc_type,scan_num)=reshape([datastruct(scan_num,proc_type).cifit_result.area_Esk_grains],syze);
            not_g(:,:,proc_type,scan_num)=reshape([datastruct(scan_num,proc_type).cifit_result.area_notindexed_grains],syze);
        
            sum_g(:,:,proc_type,scan_num)=di_g(:,:,proc_type,scan_num)+tmt_g(:,:,proc_type,scan_num)+esk_g(:,:,proc_type,scan_num)+not_g(:,:,proc_type,scan_num);
            
            % Calculate area fractions
            afrac_di_g(:,:,proc_type,scan_num)=di_g(:,:,proc_type,scan_num)./sum_g(:,:,proc_type,scan_num);
            afrac_tmt_g(:,:,proc_type,scan_num)=tmt_g(:,:,proc_type,scan_num)./sum_g(:,:,proc_type,scan_num);
            afrac_esk_g(:,:,proc_type,scan_num)=esk_g(:,:,proc_type,scan_num)./sum_g(:,:,proc_type,scan_num);
            afrac_not_g(:,:,proc_type,scan_num)=not_g(:,:,proc_type,scan_num)./sum_g(:,:,proc_type,scan_num);
        
        end
        
        % jj gives column number of the Fit or CI threshold matrix to be
        % plotted, and selects corresponding values to be plotted on y-axis
        % two lines are plotted, one for standardised, one for
        % non-standardised scan
        switch jj
            case 1
            plot(cithresh(:,jj),100*(afrac_di_g(:,jj,2,scan_num)),'marker','o','Color',color(scan_num,:),'MarkerFaceColor','w','LineStyle',':','LineWidth',lw,'MarkerSize',msz)    
            hold on
            plot(cithresh(:,jj),100*(afrac_di_g(:,jj,1,scan_num)),'marker','x','Color',color(scan_num,:),'LineWidth',lw,'MarkerSize',msz)
            xlim([0.04 0.16])
            xticks([0.05:0.025:0.15])
            xlabel('CI thresh. (only)')
            xtickformat('%-5.3f')
            ylabel(ylab)
            case 2 
            plot(fitthresh(:,jj),100*(afrac_di_g(:,jj,2,scan_num)),'marker','o','Color',color(scan_num,:),'MarkerFaceColor','w','LineStyle',':','LineWidth',lw,'MarkerSize',msz)
            hold on
            plot(fitthresh(:,jj),100*(afrac_di_g(:,jj,1,scan_num)),'marker','x','Color',color(scan_num,:),'LineWidth',lw,'MarkerSize',msz)                    
            xlim([1.1 1.7])
            xticks([1.2:0.1:1.6])
            set ( gca, 'XDir', 'reverse' )
            xlabel('Fit thresh. (only)');
            xtickformat('%-5.2f')
            ylabel(ylab)
            case 3    
            plot(cithresh(:,jj),100*(afrac_di_g(:,jj,2,scan_num)),'marker','o','Color',color(scan_num,:),'MarkerFaceColor','w','LineStyle',':','LineWidth',lw,'MarkerSize',msz)
            hold on
            plot(cithresh(:,jj),100*(afrac_di_g(:,jj,1,scan_num)),'marker','x','Color',color(scan_num,:),'LineWidth',lw,'MarkerSize',msz)      
            xlim([0.04 0.16])
            xticks([0.05:0.025:0.15])
            xlabel('CI thresh. (Fit = 1.4)')
            xtickformat('%-5.3f')
            ylabel(ylab)
            case 4    
            plot(fitthresh(:,jj),100*(afrac_di_g(:,jj,2,scan_num)),'marker','o','Color',color(scan_num,:),'MarkerFaceColor','w','LineStyle',':','LineWidth',lw,'MarkerSize',msz)
            hold on
            plot(fitthresh(:,jj),100*(afrac_di_g(:,jj,1,scan_num)),'marker','x','Color',color(scan_num,:),'LineWidth',lw,'MarkerSize',msz)                       
            xlim([1.1 1.7])
            xticks([1.2:0.1:1.6])
            set ( gca, 'XDir', 'reverse' )
            xlabel('Fit thresh. (CI = 0.1)')
            xtickformat('%-5.2f')
            ylabel(ylab)
        end
        ylim([30 70])
        
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

end

cd(savelocation)

f = gcf;
f.Position = [leftfrac*w botfrac*h horzfrac*w horzfrac*vert_div_horiz*w];

exportgraphics(f,'Cpx_crystallinity_CIfit.png','Resolution',300)

%% Fig 2b

ylab='\phi_T_m_t (%)';

% titlelist={'CI thresh only','Fit thresh only','CI thresh (Fit=1.4)','Fit Thresh (CI=0.1)'};

figure,
% For 2 subplots...
t=tiledlayout(1,2);

for jj=1:2

nexttile    
% subplot(1,4,jj)

    % Loop over each scan and processing type
    for scan_num=1:4
        for proc_type=1:2      
        
            % For each scan & proc type, obtain a 5x4 matrix of results the
            % same size as the matrix of input CI and Fit combinations,
            % store these all together in one indexable 4D matrix
            di_g(:,:,proc_type,scan_num)=reshape([datastruct(scan_num,proc_type).cifit_result.area_Di_grains],syze);
            tmt_g(:,:,proc_type,scan_num)=reshape([datastruct(scan_num,proc_type).cifit_result.area_Tmt_grains],syze);
            esk_g(:,:,proc_type,scan_num)=reshape([datastruct(scan_num,proc_type).cifit_result.area_Esk_grains],syze);
            not_g(:,:,proc_type,scan_num)=reshape([datastruct(scan_num,proc_type).cifit_result.area_notindexed_grains],syze);
        
            sum_g(:,:,proc_type,scan_num)=di_g(:,:,proc_type,scan_num)+tmt_g(:,:,proc_type,scan_num)+esk_g(:,:,proc_type,scan_num)+not_g(:,:,proc_type,scan_num);
            
            % Calculate area fractions
            afrac_di_g(:,:,proc_type,scan_num)=di_g(:,:,proc_type,scan_num)./sum_g(:,:,proc_type,scan_num);
            afrac_tmt_g(:,:,proc_type,scan_num)=tmt_g(:,:,proc_type,scan_num)./sum_g(:,:,proc_type,scan_num);
            afrac_esk_g(:,:,proc_type,scan_num)=esk_g(:,:,proc_type,scan_num)./sum_g(:,:,proc_type,scan_num);
            afrac_not_g(:,:,proc_type,scan_num)=not_g(:,:,proc_type,scan_num)./sum_g(:,:,proc_type,scan_num);
        
        end
        
        % jj gives column number of the Fit or CI threshold matrix to be
        % plotted, and selects corresponding values to be plotted on y-axis
        % two lines are plotted, one for standardised, one for
        % non-standardised scan
        switch jj
            case 1
            plot(cithresh(:,jj),100*(afrac_tmt_g(:,jj,2,scan_num)),'marker','o','Color',color(scan_num,:),'MarkerFaceColor','w','LineStyle',':','LineWidth',lw,'MarkerSize',msz)    
            hold on
            plot(cithresh(:,jj),100*(afrac_tmt_g(:,jj,1,scan_num)),'marker','x','Color',color(scan_num,:),'LineWidth',lw,'MarkerSize',msz)
            xlim([0.04 0.16])
            xticks([0.05:0.025:0.15])
            xlabel('CI thresh. (only)')
            xtickformat('%-5.3f')
            ylabel(ylab)
            case 2 
            plot(fitthresh(:,jj),100*(afrac_tmt_g(:,jj,2,scan_num)),'marker','o','Color',color(scan_num,:),'MarkerFaceColor','w','LineStyle',':','LineWidth',lw,'MarkerSize',msz)
            hold on
            plot(fitthresh(:,jj),100*(afrac_tmt_g(:,jj,1,scan_num)),'marker','x','Color',color(scan_num,:),'LineWidth',lw,'MarkerSize',msz)                    
            xlim([1.1 1.7])
            xticks([1.2:0.1:1.6])
            set ( gca, 'XDir', 'reverse' )
            xlabel('Fit thresh. (only)');
            xtickformat('%-5.2f')
            ylabel(ylab)

        end
        ylim([1 4])
        
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

end

cd(savelocation)

f = gcf;
f.Position = [leftfrac*w botfrac*h 0.98*(horzfrac/2)*w 1*horzfrac*vert_div_horiz*w];

exportgraphics(f,'Tmt_crystallinity_CIfit.png','Resolution',300)
%% Fig 2c

ylab='L_m_a_x Cpx (\mum)';

figure,
% For 2 subplots...
t=tiledlayout(1,2);

for jj=1:2

nexttile    
    for scan_num=1:4
        for proc_type=1:2      
        
            % cifit_result doesn't contain pre-calculated Lmax10, need to
            % caclulate it here...
            for ii=1:length(cithresh(:,jj))
                majoraxD=datastruct(scan_num,proc_type).cifit_result(ii,jj).majoraxD;
                [~,sort_areaID]=sort(datastruct(scan_num,proc_type).cifit_result(ii,jj).ci_fit_grains('Diopside').area);
                majoraxD=majoraxD(sort_areaID);
                lmax10_D(ii)=mean([majoraxD(end-9:end)]);
            end
            means(:,jj,proc_type,scan_num)=lmax10_D;
        
        end
        
        switch jj
            case 1    
            plot(cithresh(:,jj),means(:,jj,2,scan_num),'marker','o','Color',color(scan_num,:),'MarkerFaceColor','w','LineStyle',':','LineWidth',lw,'MarkerSize',msz)
            hold on
            plot(cithresh(:,jj),means(:,jj,1,scan_num),'marker','x','Color',color(scan_num,:),'LineWidth',lw,'MarkerSize',msz)                     
                    xlim([0.04 0.16])
                    xticks([0.05:0.025:0.15])
            xlabel('CI thresh. (only)')
            xtickformat('%-5.3f')
            ylabel(ylab)
            case 2 
            plot(fitthresh(:,jj),means(:,jj,2,scan_num),'marker','o','Color',color(scan_num,:),'MarkerFaceColor','w','LineStyle',':','LineWidth',lw,'MarkerSize',msz)
            hold on
            plot(fitthresh(:,jj),means(:,jj,1,scan_num),'marker','x','Color',color(scan_num,:),'LineWidth',lw,'MarkerSize',msz)
                    xlim([1.1 1.7])
                    xticks([1.2:0.1:1.6])
            set ( gca, 'XDir', 'reverse' )
            fitlabel=xlabel('Fit thresh. (only)');
            xtickformat('%-5.2f')
            ylabel(ylab)
        end

        ylim([50 125])
    
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

end

cd(savelocation)

f = gcf;
f.Position = [leftfrac*w botfrac*h 1*(horzfrac/2)*w 1*horzfrac*vert_div_horiz*w];

exportgraphics(f,'Cpx_Lmax_CIfit.png','Resolution',300)
%% Fig 2d

ylab='S_v^P Cpx (mm^-^1)';

figure,
% For 2 subplots...
t=tiledlayout(1,2);

for jj=1:2

nexttile    
    for scan_num=1:4
        for proc_type=1:2      
        
                gbldi(:,:,proc_type,scan_num)=reshape([datastruct(scan_num,proc_type).cifit_result.GBlength_Di],syze);
            
                di_g(:,:,proc_type,scan_num)=reshape([datastruct(scan_num,proc_type).cifit_result.area_Di_grains],syze);
                tmt_g(:,:,proc_type,scan_num)=reshape([datastruct(scan_num,proc_type).cifit_result.area_Tmt_grains],syze);
                esk_g(:,:,proc_type,scan_num)=reshape([datastruct(scan_num,proc_type).cifit_result.area_Esk_grains],syze);
                not_g(:,:,proc_type,scan_num)=reshape([datastruct(scan_num,proc_type).cifit_result.area_notindexed_grains],syze);
            
                sum_g(:,:,proc_type,scan_num)=di_g(:,:,proc_type,scan_num)+tmt_g(:,:,proc_type,scan_num)+esk_g(:,:,proc_type,scan_num)+not_g(:,:,proc_type,scan_num);
                
                % %
                afrac_di_g(:,:,proc_type,scan_num)=di_g(:,:,proc_type,scan_num)./sum_g(:,:,proc_type,scan_num);
                afrac_tmt_g(:,:,proc_type,scan_num)=tmt_g(:,:,proc_type,scan_num)./sum_g(:,:,proc_type,scan_num);
                afrac_esk_g(:,:,proc_type,scan_num)=esk_g(:,:,proc_type,scan_num)./sum_g(:,:,proc_type,scan_num);
                afrac_not_g=not_g(:,:,proc_type,scan_num)./sum_g(:,:,proc_type,scan_num);
            
                svp(:,:,proc_type,scan_num)=((4/pi())*(gbldi(:,:,proc_type,scan_num)/1000)./...
                (datastruct(scan_num).true_area/1e6))./afrac_di_g(:,:,proc_type,scan_num);  
        
        end
        
        switch jj
            case 1    
            plot(cithresh(:,jj),svp(:,jj,2,scan_num),'marker','o','Color',color(scan_num,:),'MarkerFaceColor','w','LineStyle',':','LineWidth',lw,'MarkerSize',msz)
            hold on
            plot(cithresh(:,jj),svp(:,jj,1,scan_num),'marker','x','Color',color(scan_num,:),'LineWidth',lw,'MarkerSize',msz)                     
                    xlim([0.04 0.16])
                    xticks([0.05:0.025:0.15])
            xlabel('CI thresh. (only)')
            xtickformat('%-5.3f')
            ylabel(ylab)
            case 2 
            plot(fitthresh(:,jj),svp(:,jj,2,scan_num),'marker','o','Color',color(scan_num,:),'MarkerFaceColor','w','LineStyle',':','LineWidth',lw,'MarkerSize',msz)
            hold on
            plot(fitthresh(:,jj),svp(:,jj,1,scan_num),'marker','x','Color',color(scan_num,:),'LineWidth',lw,'MarkerSize',msz)
                    xlim([1.1 1.7])
                    xticks([1.2:0.1:1.6])
            set ( gca, 'XDir', 'reverse' )
            fitlabel=xlabel('Fit thresh. (only)');
            xtickformat('%-5.2f')
            ylabel(ylab)
        end

        ylim([200 1700])
    
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

end

cd(savelocation)

f = gcf;
f.Position = [leftfrac*w botfrac*h 1*(horzfrac/2)*w 1*horzfrac*vert_div_horiz*w];

exportgraphics(f,'Cpx_SvP_CIfit.png','Resolution',300)
%% Fig 2e

ylab='N_A Tmt (mm^-^2)';

figure,
% For 2 subplots...
t=tiledlayout(1,2);

for jj=1:2

nexttile    
    for scan_num=1:4
        for proc_type=1:2      
        
            dens(:,:,proc_type,scan_num)=reshape([datastruct(scan_num,proc_type).cifit_result.N_density_Mt_true],syze);
        
        end
        
        switch jj
            case 1    
            plot(cithresh(:,jj),dens(:,jj,2,scan_num),'marker','o','Color',color(scan_num,:),'MarkerFaceColor','w','LineStyle',':','LineWidth',lw,'MarkerSize',msz)
            hold on
            plot(cithresh(:,jj),dens(:,jj,1,scan_num),'marker','x','Color',color(scan_num,:),'LineWidth',lw,'MarkerSize',msz)                     
                    xlim([0.04 0.16])
                    xticks([0.05:0.025:0.15])
            xlabel('CI thresh. (only)')
            xtickformat('%-5.3f')
            ylabel(ylab)
            case 2 
            plot(fitthresh(:,jj),dens(:,jj,2,scan_num),'marker','o','Color',color(scan_num,:),'MarkerFaceColor','w','LineStyle',':','LineWidth',lw,'MarkerSize',msz)
            hold on
            plot(fitthresh(:,jj),dens(:,jj,1,scan_num),'marker','x','Color',color(scan_num,:),'LineWidth',lw,'MarkerSize',msz)
                    xlim([1.1 1.7])
                    xticks([1.2:0.1:1.6])
            set ( gca, 'XDir', 'reverse' )
            fitlabel=xlabel('Fit thresh. (only)');
            xtickformat('%-5.2f')
            ylabel(ylab)
        end

        ylim([1000 7000])
    
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

end

cd(savelocation)

f = gcf;
f.Position = [leftfrac*w botfrac*h 1*(horzfrac/2)*w 1*horzfrac*vert_div_horiz*w];

exportgraphics(f,'Tmt_NA_CIfit.png','Resolution',300)
