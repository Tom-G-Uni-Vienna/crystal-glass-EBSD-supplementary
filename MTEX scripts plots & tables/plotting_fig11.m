%% Load data & setup formatting

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
input1=load('ETK5a5b_allscans_distort_regrid_final');

datastruct=input1.datastruct;

clear input1 

fprintf('MATLAB variables loaded!\n\n') 

%%
% Cpx crystal symmetry as used for indexing
CSdi=crystalSymmetry('2/m', [9.585 8.776 5.26], [90,106.85,90]*degree,...
    'X||a*', 'Y||b*', 'Z||c', 'mineral', 'Diopside', 'color', [0.8 0.8 0.8]);

% Helper variable for the combination of data - needed as 
% sample 1 = maps 1+2, sample 2 = maps 3+4 ->
% so the two maps for a given sample in the input variable datastruct are
% at row positions q(sample_num) and q(sample_num)+1.
q=[1,3]; 

combined_cpx_oris=cell(2,1);
odf20=cell(2,1);

for sample_num=1:2
    proc_type=2;
        combined_cpx_oris{sample_num,1} =...
            [datastruct(q(sample_num),proc_type).grains2('Diopside').meanOrientation;...
            datastruct(q(sample_num)+1,proc_type).grains2('Diopside').meanOrientation ];

        odf20{sample_num,1}=calcDensity(combined_cpx_oris{sample_num,1},'halfwidth',20*degree);

    
end
%%
pfAnnotations = @(varargin) [];
setMTEXpref('pfAnnotations',pfAnnotations);

screendat=get(0,'ScreenSize');
w=screendat(3); % width
h=screendat(4); % height
leftfrac=0.02;
botfrac=0.02;
horzfrac=0.4;
vertfrac=0.8;

figure
fig2=mtexFigure('layout',size(combined_cpx_oris));
for sample_num=1:2
     
    nextAxis
    plotPDF(odf20{sample_num},Miller(0,1,0,CSdi,'hkl'),...
            'projection','eangle','upper','antipodal')

    textureindex(odf20{sample_num})
   

end
mtexColorbar
colormap("parula")


cd(savelocation)

f = gcf;
f.Position = [leftfrac*w botfrac*h horzfrac*w vertfrac*h];
exportgraphics(f,'odfs_cpx010.png','Resolution',400)


