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

fprintf('loading MATLAB variables...\n\n') 

cd(loadlocation)

% MATLAB prefers to load data to structures and extract variables from them
% instead of loading the variables directly 
% (https://de.mathworks.com/matlabcentral/answers/28676-why-use-x-load-myfile-mat)
input1=load('ETK5a5b_allscans_distort.mat','datastruct');

datastruct=input1.datastruct;

clear input1

fprintf('MATLAB variables loaded!\n\n') 

%%

screendat=get(0,'ScreenSize');
w=screendat(3);
h=screendat(4);

leftfrac=0.02;
botfrac=0.045;
horzfrac=0.55;

ar = 0.9;

lw=2.5; % linewidth
fs=28; % fontsize x and y label
ax_fs = 23; % fontsize axis numbering
ax_lw = 2.5; % axis linewidth

scanname={'Anhy\_150','Anhy\_200','Hy\_200','Hy\_1000'};

figure,
tiledlayout(2,2)
for i=1:4
    nexttile
    fits=datastruct(i,1).data('indexed').fit;
    histogram(fits(datastruct(i,1).data('indexed').ci>0.1),100)
    xlim([0.4 2])
    xticks([0.4 1 1.4 2])
    title(scanname{i})
    ylabel('N pixels')
    xlabel('Fit')
        ax = gca;
    ax.XAxis.FontSize = ax_fs;
    ax.YAxis.FontSize = ax_fs;
    ax.XLabel.FontSize = fs-3;
    ax.YLabel.FontSize = fs;
    ax.Title.FontSize = fs;
    ax.LineWidth = ax_lw;

end

cd(savelocation)

f = gcf;
f.Position = [leftfrac*w botfrac*h horzfrac*w horzfrac*w*ar];
exportgraphics(f,'fit_histograms.png','Resolution',450)
