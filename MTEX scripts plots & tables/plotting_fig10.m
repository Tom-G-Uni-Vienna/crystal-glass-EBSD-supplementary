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
    '\FINAL Supplementary\FINAL scripts & data\MATLAB data\Data undistorted maps'];

% Output path (automatically derived as long as supp. mat. structure
% remains unchanged)
savelocation = [ supplementary_location...
    '\FINAL Supplementary\FINAL scripts & data\Figures & tables' ];

fprintf('loading MATLAB variables...\n\n') 

cd(loadlocation)

% MATLAB prefers to load data to structures and extract variables from them
% instead of loading the variables directly 
% (https://de.mathworks.com/matlabcentral/answers/28676-why-use-x-load-myfile-mat)
input1=load('ETK5a5b_allscans_NOdistort_cithresh.mat','datastruct');

datastruct=input1.datastruct;

clear input1

fprintf('MATLAB variables loaded!\n\n') 

% setMTEXpref('showMicronBar','off')

%%
setMTEXpref('showCoordinates','on')

screendat=get(0,'ScreenSize');
w=screendat(3);
h=screendat(4);

leftfrac=0.02;
botfrac=0.045;
horzfrac=0.32;

rect=[45.9 -218.5 9.4 9.4];

cd(savelocation)

lw=3;

cithresh=datastruct(1).cithresh;

mapn=2;

ar=rect(4)/rect(3);

cifit_col=1;

color_N=256;

% color = colororder('glow')
color = hot(color_N);
% color = jet(color_N);
% color = copper(color_N);
% color = pink(color_N);
color = flipud(color);
color=color(round(linspace(1,color_N,length(cithresh))),:);

ebsd=datastruct(mapn,1).data;

sel=ebsd.inpolygon(rect);

ebsdsel=ebsd(sel).gridify;

figure,
plot(ebsdsel,ebsdsel.iq)
mtexColorMap black2white
clim([20000 120000])
hold on
for i = 1:5
    hold on
  % s=smooth(datastruct(mapn).cifit_result(i,cifit_col).ci_fit_grains,15);  
  s=datastruct(mapn).cifit_result(i,cifit_col).ci_fit_grains;
  plot(s.boundary,'linewidth',lw,'linecolor',color(i,:))
  clear s
end
%     hold on
%   plot(datastruct(mapn).cifit_result(2,cifit_col).ci_fit_grains('Diopside').boundary,'linewidth',1.5,'linecolor','b')
%     hold on
%   plot(datastruct(mapn).cifit_result(4,cifit_col).ci_fit_grains('Diopside').boundary,'linewidth',1.5,'linecolor','r')
%     hold on
%   plot(datastruct(mapn).cifit_result(3,cifit_col).ci_fit_grains('Diopside').boundary,'linewidth',1.5,'linecolor','k')
% hold off
xlim([rect(1) rect(1)+rect(3)])
ylim([rect(2) rect(2)+rect(4)])

f = gcf;
f.Position = [leftfrac*w botfrac*h horzfrac*w horzfrac*w*ar];
% exportgraphics(f,'stand_v_nostand_IQ_hotcolors_smooth.png','Resolution',300)
exportgraphics(f,'stand_v_nostand_IQ_hotcolors.png','Resolution',300)

horzfrac=0.9

% Area grains-based HOT colors
syze=size(datastruct(1).cifit_result);

lw=2.5
fs=23
afs=21
ms=15

figure,
for jj=3

subplot(1,4,jj)

    for ee=mapn
        for ff=1:2      
        
            di_g(:,:,ff,ee)=reshape([datastruct(ee,ff).cifit_result.area_Di_grains],syze);
            tmt_g(:,:,ff,ee)=reshape([datastruct(ee,ff).cifit_result.area_Tmt_grains],syze);
            esk_g(:,:,ff,ee)=reshape([datastruct(ee,ff).cifit_result.area_Esk_grains],syze);
            not_g(:,:,ff,ee)=reshape([datastruct(ee,ff).cifit_result.area_notindexed_grains],syze);
        
            sum_g(:,:,ff,ee)=di_g(:,:,ff,ee)+tmt_g(:,:,ff,ee)+esk_g(:,:,ff,ee)+not_g(:,:,ff,ee);
            
            % %
            afrac_di_g(:,:,ff,ee)=di_g(:,:,ff,ee)./sum_g(:,:,ff,ee);
            afrac_tmt_g(:,:,ff,ee)=tmt_g(:,:,ff,ee)./sum_g(:,:,ff,ee);
            afrac_esk_g(:,:,ff,ee)=esk_g(:,:,ff,ee)./sum_g(:,:,ff,ee);
            afrac_not_g(:,:,ff,ee)=not_g(:,:,ff,ee)./sum_g(:,:,ff,ee);
        
        end

        plot(cithresh(:,jj),100*(1-afrac_not_g(:,jj,1,ee)),'marker','o','Color','k','LineWidth',lw,'MarkerEdgeColor','none','MarkerFaceColor','none','MarkerSize',ms)
        hold on
        for qq=1:5

        plot(cithresh(qq,jj),100*(1-afrac_not_g(qq,jj,1,ee)),'marker','o','Color','none','LineWidth',lw,'MarkerEdgeColor','k','MarkerFaceColor',color(qq,:),'MarkerSize',ms)
        hold on
        end     
            xlim([0.04 0.16])
            xticks([0.05:0.025:0.15])
            xlabel('CI Thresh (Fit = 1.4)')

        ylim([57 62])
        
        ylabel('\phi_c_r_y_s_t_a_l_s (%)')

        ax = gca
        ax.XAxis.FontSize = afs;
        ax.YAxis.FontSize = afs;
        ax.XLabel.FontSize = fs;
        ax.YLabel.FontSize = fs;
        ax.LineWidth = 2.5
    end
end

f = gcf;
f.Position = [leftfrac*w botfrac*h horzfrac*w horzfrac*w*0.4];
exportgraphics(f,'crystallinity_hotcolors.png','Resolution',300)

%%
setMTEXpref('showCoordinates','on')

screendat=get(0,'ScreenSize');
w=screendat(3);
h=screendat(4);

leftfrac=0.02;
botfrac=0.045;
horzfrac=0.32;

rect=[189 -222 9 9];

cd(savelocation)

lw=3;

cithresh=datastruct(1).cithresh;

mapn=1;

ar=rect(4)/rect(3);

cifit_col=1;

color_N=256;

% color = colororder('glow')
color = hot(color_N);
% color = jet(color_N);
% color = copper(color_N);
% color = pink(color_N);
color = flipud(color);
color=color(round(linspace(1,color_N,length(cithresh))),:);

ebsd=datastruct(mapn,1).data;

sel=ebsd.inpolygon(rect);

ebsdsel=ebsd(sel).gridify;

figure,
plot(ebsdsel,ebsdsel.iq)
mtexColorMap black2white
clim([20000 120000])
hold on
for i = 1:5
    hold on
  % s=smooth(datastruct(mapn).cifit_result(i,cifit_col).ci_fit_grains,15);  
  s=datastruct(mapn).cifit_result(i,cifit_col).ci_fit_grains;
  plot(s.boundary,'linewidth',lw,'linecolor',color(i,:))
  clear s
end
%     hold on
%   plot(datastruct(mapn).cifit_result(2,cifit_col).ci_fit_grains('Diopside').boundary,'linewidth',1.5,'linecolor','b')
%     hold on
%   plot(datastruct(mapn).cifit_result(4,cifit_col).ci_fit_grains('Diopside').boundary,'linewidth',1.5,'linecolor','r')
%     hold on
%   plot(datastruct(mapn).cifit_result(3,cifit_col).ci_fit_grains('Diopside').boundary,'linewidth',1.5,'linecolor','k')
% hold off
xlim([rect(1) rect(1)+rect(3)])
ylim([rect(2) rect(2)+rect(4)])

f = gcf;
f.Position = [leftfrac*w botfrac*h horzfrac*w horzfrac*w*ar];
% exportgraphics(f,'stand_v_nostand_IQ_hotcolors_smooth.png','Resolution',300)
exportgraphics(f,'stand_v_nostand_IQ_hotcolors2.png','Resolution',300)

horzfrac=0.9

% Area grains-based HOT colors
syze=size(datastruct(1).cifit_result);

lw=2.5
fs=23
afs=21
ms=15

figure,
for jj=3

subplot(1,4,jj)

    for ee=mapn
        for ff=1:2      
        
            di_g(:,:,ff,ee)=reshape([datastruct(ee,ff).cifit_result.area_Di_grains],syze);
            tmt_g(:,:,ff,ee)=reshape([datastruct(ee,ff).cifit_result.area_Tmt_grains],syze);
            esk_g(:,:,ff,ee)=reshape([datastruct(ee,ff).cifit_result.area_Esk_grains],syze);
            not_g(:,:,ff,ee)=reshape([datastruct(ee,ff).cifit_result.area_notindexed_grains],syze);
        
            sum_g(:,:,ff,ee)=di_g(:,:,ff,ee)+tmt_g(:,:,ff,ee)+esk_g(:,:,ff,ee)+not_g(:,:,ff,ee);
            
            % %
            afrac_di_g(:,:,ff,ee)=di_g(:,:,ff,ee)./sum_g(:,:,ff,ee);
            afrac_tmt_g(:,:,ff,ee)=tmt_g(:,:,ff,ee)./sum_g(:,:,ff,ee);
            afrac_esk_g(:,:,ff,ee)=esk_g(:,:,ff,ee)./sum_g(:,:,ff,ee);
            afrac_not_g(:,:,ff,ee)=not_g(:,:,ff,ee)./sum_g(:,:,ff,ee);
        
        end

        plot(cithresh(:,jj),100*(1-afrac_not_g(:,jj,1,ee)),'marker','o','Color','k','LineWidth',lw,'MarkerEdgeColor','none','MarkerFaceColor','none','MarkerSize',ms)
        hold on
        for qq=1:5

        plot(cithresh(qq,jj),100*(1-afrac_not_g(qq,jj,1,ee)),'marker','o','Color','none','LineWidth',lw,'MarkerEdgeColor','k','MarkerFaceColor',color(qq,:),'MarkerSize',ms)
        hold on
        end     
            xlim([0.04 0.16])
            xticks([0.05:0.025:0.15])
            xlabel('CI Thresh (Fit = 1.4)')

        ylim([47 55])
        
        ylabel('\phi_c_r_y_s_t_a_l_s (%)')

        ax = gca
        ax.XAxis.FontSize = afs;
        ax.YAxis.FontSize = afs;
        ax.XLabel.FontSize = fs;
        ax.YLabel.FontSize = fs;
        ax.LineWidth = 2.5
    end
end

f = gcf;
f.Position = [leftfrac*w botfrac*h horzfrac*w horzfrac*w*0.4];
exportgraphics(f,'crystallinity_hotcolors2.png','Resolution',300)
