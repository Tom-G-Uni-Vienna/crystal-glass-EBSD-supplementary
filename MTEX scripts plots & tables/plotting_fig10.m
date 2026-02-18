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

lw=3; % linewidth
fs=21; % fontsize x and y label

% Get screen sizes and define screen fractions for plot positioning
screendat=get(0,'ScreenSize');
w=screendat(3); % width
h=screendat(4); % height
leftfrac=0.02;
botfrac=0.02;
horzfrac=0.4;
vertfrac=0.85;

% Generate color list of desired order
color=colororder('gem');
c2=color(2,:);
c4=color(4,:);
color(4,:)=c2;
color(2,:)=c4;
color=flip(color);
color=color([4:7],:);


%% Rose diagrams of boundary segment directions

% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

% Choose a particular scan for which to plot rose diagrams
scan_num=1;
proc_type=1;

% Update list of smoothing numbers to add 0 at the beginning
smoothnum2=[0,smoothnum];

% Create variable containing boundary vectors for unsmoothed PLUS all
% smoothing iteratations
gbdirs_plus_unsmooth=cell(length(smoothnum2),1);
gbdirs_plus_unsmooth{1}=datastruct(scan_num,proc_type).grains2.boundary.direction;
for i = 1:12
  gbdirs_plus_unsmooth{i+1} = datastruct(scan_num).smooth_result(i).smooth_grains.boundary.direction;
end

% Sizes of rose diagrams 
rlim_list=[0.18,0.07,0.07,0.07,0.07,0.07];

% Plot frequency normalised rose diagrams of boundary directions for some,
% but not all of the N smoothing iterations tested, with only the first
% diagram having a different scaling to the rest
figure
inc = 0;
for i=[0,5,6,10,11,12]
    inc=inc+1;
  subplot(3,2,inc)
  % 'weights', norm() means that rose diagram is NOT length weighted (norms
  % of all vectors are 1)
  histogram(gbdirs_plus_unsmooth{i+1}, 'weights',norm(gbdirs_plus_unsmooth{i+1}),30,'EdgeColor','none','FaceColor','b','FaceAlpha',1)
  title([ num2str(smoothnum2(i+1)) ' iterations'])
  ax=gca;
    ax.RTick = [0.07,0.18];
    ax.RTickLabel = {};
    thetaticks([0 90 180 270]);
    ax.RLim = [0 rlim_list(inc)];
    ax.LineWidth = lw;
    ax.Title.FontSize = fs;

end

cd(savelocation)

f = gcf;
f.Position = [leftfrac*w botfrac*h horzfrac*w vertfrac*w];
exportgraphics(f,'smoothnum_rosediagram.png','Resolution',400)
% close


