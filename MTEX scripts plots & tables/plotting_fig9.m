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

mapn=2;

ar=rect(4)/rect(3);

ebsd=datastruct(mapn,1).data;

sel=ebsd.inpolygon(rect);

ebsdsel=ebsd(sel).gridify;

figure,
plot(ebsdsel,ebsdsel.iq)
mtexColorMap black2white
clim([20000 120000])
hold on
plot(datastruct(mapn).cifit_result(3,1).ci_fit_grains.boundary,'linewidth',lw,'linecolor','HotPink')
hold on
plot(datastruct(mapn+4).cifit_result(3,1).ci_fit_grains.boundary,'linewidth',lw,'linecolor','Orange')
xlim([rect(1) rect(1)+rect(3)])
ylim([rect(2) rect(2)+rect(4)])

f = gcf;
f.Position = [leftfrac*w botfrac*h horzfrac*w horzfrac*w*ar];
exportgraphics(f,'stand_v_nostand_gbs_IQ.png','Resolution',300)

figure,
plot(ebsdsel,ebsdsel.ci)
mtexColorMap black2white
clim([0 1])
hold on
plot(datastruct(mapn).cifit_result(3,1).ci_fit_grains.boundary,'linewidth',lw,'linecolor','HotPink')
hold on
plot(datastruct(mapn+4).cifit_result(3,1).ci_fit_grains.boundary,'linewidth',lw,'linecolor','Orange')
xlim([rect(1) rect(1)+rect(3)])
ylim([rect(2) rect(2)+rect(4)])

f = gcf;
f.Position = [leftfrac*w botfrac*h horzfrac*w horzfrac*w*ar];
exportgraphics(f,'stand_v_nostand_gbs_nostandCI.png','Resolution',300)

ebsd=datastruct(mapn,2).data;

sel=ebsd.inpolygon(rect);

ebsdsel=ebsd(sel).gridify;

figure,
plot(ebsdsel,ebsdsel.ci)
mtexColorMap black2white
clim([0 1])
hold on
plot(datastruct(mapn).cifit_result(3,1).ci_fit_grains.boundary,'linewidth',lw,'linecolor','HotPink')
hold on
plot(datastruct(mapn+4).cifit_result(3,1).ci_fit_grains.boundary,'linewidth',lw,'linecolor','Orange')
xlim([rect(1) rect(1)+rect(3)])
ylim([rect(2) rect(2)+rect(4)])

f = gcf;
f.Position = [leftfrac*w botfrac*h horzfrac*w horzfrac*w*ar];
exportgraphics(f,'stand_v_nostand_gbs_standCI.png','Resolution',300)

ebsdsel_di=ebsdsel('Diopside');

ipf_sel=ipfHSVKey(ebsdsel_di.orientations)
ipf_sel.inversePoleFigureDirection=zvector;
ipf_sel_colors=ipf_sel.orientation2color(ebsdsel_di.orientations)


figure,
plot(ebsdsel,ebsdsel.ci)
mtexColorMap black2white
clim([0 1])
hold on
plot(ebsdsel_di(ebsdsel_di.ci>0.1),ipf_sel_colors(ebsdsel_di.ci>0.1,:))
hold on
plot(datastruct(mapn).cifit_result(3,1).ci_fit_grains.boundary,'linewidth',lw,'linecolor','HotPink')
hold on
plot(datastruct(mapn+4).cifit_result(3,1).ci_fit_grains.boundary,'linewidth',lw,'linecolor','Orange')
xlim([rect(1) rect(1)+rect(3)])
ylim([rect(2) rect(2)+rect(4)])

f = gcf;
f.Position = [leftfrac*w botfrac*h horzfrac*w horzfrac*w*ar];
exportgraphics(f,'stand_v_nostand_gbs_standIPF.png','Resolution',300)




