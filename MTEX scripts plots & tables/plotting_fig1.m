%% Load data & setup formatting

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
input1=load('ETK5a5b_allscans_NOdistort.mat','datastruct');

datastruct=input1.datastruct;

clear input1 

fprintf('MATLAB variables loaded!\n\n') 

%%
setMTEXpref('showCoordinates','on')
setMTEXpref('showMicronBar','off')

screendat=get(0,'ScreenSize');
w=screendat(3);
h=screendat(4);
aratio=1;

leftfrac=0.02;
botfrac=0.045;
horzfrac=0.6;

for mapn=1:4

ebsd=datastruct(mapn,2).data;

ebsd(ebsd.ci<0.1).phase=0;
ebsd(ebsd.fit>1.45).phase=0;

figure,plot(ebsd,ebsd.iq)
mtexColorMap black2white
hold on
plot(ebsd('Magnetite'),'FaceAlpha',0.2,'FaceColor','r')
if ~isempty(ebsd('Eskolaite'))
hold on
plot(ebsd('Eskolaite'),'FaceAlpha',0.2,'FaceColor','LightBlue')
end
ax=gca;
legend(ax,'off')

cd(savelocation)

% f = gcf;
% f.Position = [leftfrac*w botfrac*h horzfrac*w horzfrac*w*aratio];
% exportgraphics(f,[ datastruct(mapn,2).mapname(2:end-9) '_IQ_undistort.png'],'Resolution',300)

clear ebsd

% close(f)

end
