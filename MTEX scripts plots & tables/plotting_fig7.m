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
input1=load('ETK5a5b_allscans_distort_regrid_final.mat','datastruct');
input2=load('ETK5a5b_allscans_distort_combined_result.mat','resultstruct');

d=input1.datastruct;
r=input2.resultstruct;

clear input1 input2

fprintf('MATLAB variables loaded!\n\n') 

% Helper variable for the combination of data - needed as 
% sample 1 = maps 1+2, sample 2 = maps 3+4 ->
% so the two maps for a given sample in the input variable datastruct are
% at row positions q(sample_num) and q(sample_num)+1.
q=[1,3]; 

%% Pontesilli et al. 2019 data (10.1016/j.chemgeo.2019.02.015)
% Data taken from tables 2 & 3 (with corrections as mentioned in
% supplementary text of this work)

% "Dry" = ETK5a, "wet" = ETK5b. 

% Diopside
% Rows with numbers correspond to following:
% Area frac, sqrt(L*W)max, Gmax, SvP. 
% (2nd column SV = sample variation or s.d. = standard deviation):
% SV, s.d., s.d., SV.
% units
% - , µm, cms^-1, mm^-1
TPDdry =  [ [0.5353  ; 0 ; 53.66  ; 0 ; 14.91e-7 ; 0    ; 0     ; 0   ; 513.84  ], ...
            [0.1241  ; 0 ; 5.39   ; 0 ; 1.50e-7 ; 0    ; 0     ; 0    ; 87.93   ]];

TPDwet =  [ [0.3846  ; 0 ; 63.24  ; 0 ; 18.02e-7 ; 0    ; 0     ; 0   ; 398.81  ], ...
            [0.0609  ; 0 ; 6.41   ; 0 ; 2.49e-7 ; 0    ; 0     ; 0    ; 92.31   ]];

% Titanomagnetite
% Rows with numbers correspond to following:
% Area frac, sqrt(L*W)max, Gmax, Gbatch, d, NA. 
% (2nd column SV = sample variation or s.d. = standard deviation):
% SV, s.d., s.d., -, -, -.
% units
% - , µm, cms^-1, cms^-1, µm, mm^-2
TPMdry =  [ [0.0077  ; 0 ; 5.31   ; 0 ; 1.48e-7 ; 0.71e-7 ; 2.56  ; 1174  ], ...
            [0.0068  ; 0 ; 0.78   ; 0 ; 0.22e-7 ; 0    ; 0     ; 0     ]];

TPMwet =  [ [0.018  ; 0 ; 7.20  ; 0 ; 2.00e-7 ; 1.13e-7 ; 4.08  ; 1083  ], ...
            [0.0157 ; 0 ; 0.56  ; 0 ; 0.16e-7 ; 0    ; 0     ; 0     ]];

TPM={TPMdry,TPMwet}
TPD={TPDdry,TPDwet}

%% Variables for figure formatting
screendat=get(0,'ScreenSize');
w=screendat(3); % screen width
h=screendat(4); % screen height

% scaling fractions of screen for fig positioning
leftfrac=0.45; 
botfrac=0.02;
horzfrac=0.33;
vertfrac=0.9;

% Control positioning of ticks along x-axis:
xmin=0;
xmax=3.25;
% x position within scatter plot for the Pont. et al. points
pos=[0.75, 2]; 
% x-shift of EBSD result point with respect to first point
posshift=0.5; 
% x shift of the individual scan results for EBSD w.r.t.
% the central "combined" result).
smallshift=0.075;

colorlist={'r','b'}; % marker colors ("dry" and "wet")

ms1=10; % markersize of "X" symbol
ms2=6; % markersize of "O" symbol
msS=8; % markersize of square symbol
lw=1; % linewidth
ax_lw = 1; % axis linewidth
axx_fs = 10; % x axis ticklabel fontsize
axy_fs = 10; % y axis ticklabel fontsize

names = {' '; 'BSE_A_n_h_y'; '  EBSD_A_n_h_y'; 'BSE_H_y'; '  EBSD_H_y'};

close
%%
figure

% Cpx SvP % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    procA=2;
    procB=1;
    row=9;
    
    if procA==1
        procMark='x';
        procMarkB='o';
        msA=ms1;
        msB=ms2;
    else
        procMark='o';
        procMarkB='x';
        msA=ms2;
        msB=ms1;
    end
    
    inputvar='SVPcpx';
    
    subplot(2,2,1)
    for sample_num=1:2
        
        sel=TPD{sample_num};

        % Pontesilli et al. data
        errorbar( pos(sample_num), sel(row,1), sel(row,2), sel(row,2) ,'Marker','s' ,'Color',colorlist{sample_num},'MarkerEdgeColor',colorlist{sample_num},'MarkerFaceColor',colorlist{sample_num},'markerSize',msS,'LineWidth',lw)
        hold on
        % EBSD data, preferred processing type (errorbar smaller than
        % symbol)
        plot( pos(sample_num)+posshift, r(sample_num,procA).(inputvar)(1), 'Marker',procMark ,'Color',colorlist{sample_num},'markerSize',msA,'LineWidth',2*lw);
        hold on
        % EBSD values for individual scans (for "preferred processing type" only)
        plot( pos(sample_num)+posshift-smallshift, d(q(sample_num),procA).(inputvar), 'Marker',procMark,'Color',colorlist{sample_num},'markerSize',msA)
        hold on
        plot( pos(sample_num)+posshift+smallshift, d(q(sample_num)+1,procA).(inputvar), 'Marker',procMark,'Color',colorlist{sample_num},'markerSize',msA)
        hold on
        % EBSD data, non-preferred processing type 
        plot( pos(sample_num)+posshift, r(sample_num,procB).(inputvar)(1), 'Marker',procMarkB,'Color',colorlist{sample_num},'markerSize',msB,'LineWidth',2*lw)
        hold on
    end
    hold off
    xlim([xmin xmax])
    ylim([250 1150])
    ylabel('Cpx S_v^P mm^-^1')
    set(gca,'xtick',[0 pos(1),pos(1)+posshift,pos(2),pos(2)+posshift],'xticklabel',names)
    ax=gca;
    ax.XTickLabelRotation = 90;
    ax.LineWidth = ax_lw;     
    ax.XAxis.FontSize = axx_fs;
    ax.YAxis.FontSize = axy_fs;

% Tmt NA % % % % % % % % % % % % % % % % % % % % % % % % % % % 

    procA=2;
    procB=1;
    row=8;
    
    if procA==1
        procMark='x';
        procMarkB='o';
        msA=ms1;
        msB=ms2;
    else
        procMark='o';
        procMarkB='x';
        msA=ms2;
        msB=ms1;
    end
    
    inputvar='NAtmt';
    
    subplot(2,2,2)
    for sample_num=1:2
        
        sel=TPM{sample_num};
        
        % Pontesilli et al. data
        errorbar( pos(sample_num), sel(row,1), sel(row,2), sel(row,2) ,'Marker','s' ,'Color',colorlist{sample_num},'MarkerEdgeColor',colorlist{sample_num},'MarkerFaceColor',colorlist{sample_num},'markerSize',msS,'LineWidth',lw)
        hold on
        % EBSD data, preferred processing type (errorbar smaller than
        % symbol)
        plot( pos(sample_num)+posshift, r(sample_num,procA).(inputvar)(1), 'Marker',procMark ,'Color',colorlist{sample_num},'markerSize',msA,'LineWidth',2*lw);
        hold on
        % EBSD values for individual scans (for "preferred processing type" only)
        plot( pos(sample_num)+posshift-smallshift, d(q(sample_num),procA).(inputvar), 'Marker',procMark,'Color',colorlist{sample_num},'markerSize',msA)
        hold on
        plot( pos(sample_num)+posshift+smallshift, d(q(sample_num)+1,procA).(inputvar), 'Marker',procMark,'Color',colorlist{sample_num},'markerSize',msA)
        hold on
        % EBSD data, non-preferred processing type 
        plot( pos(sample_num)+posshift, r(sample_num,procB).(inputvar)(1), 'Marker',procMarkB,'Color',colorlist{sample_num},'markerSize',msB,'LineWidth',2*lw)
        hold on
    end
    hold off
    xlim([xmin xmax])
    % ylim([0 4800])
    ylabel('N_A Tmt mm^-^2')
    set(gca,'xtick',[0 pos(1),pos(1)+posshift,pos(2),pos(2)+posshift],'xticklabel',names)
    ax=gca;
    ax.XTickLabelRotation = 90;
    ax.LineWidth = ax_lw;     
    ax.XAxis.FontSize = axx_fs;
    ax.YAxis.FontSize = axy_fs;

% Tmt d % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    
    procA=1;
    procB=2;
    row=7;
    
    if procA==1
        procMark='x';
        procMarkB='o';
        msA=ms1;
        msB=ms2;
    else
        procMark='o';
        procMarkB='x';
        msA=ms2;
        msB=ms1;
    end
    
    inputvar='d_Tmt';
    
    subplot(2,2,3)
    for sample_num=1:2
        
        sel=TPM{sample_num};
        
        % Pontesilli et al. data
        errorbar( pos(sample_num), sel(row,1), sel(row,2), sel(row,2) ,'Marker','s' ,'Color',colorlist{sample_num},'MarkerEdgeColor',colorlist{sample_num},'MarkerFaceColor',colorlist{sample_num},'markerSize',msS,'LineWidth',lw)
        hold on
        % EBSD data, preferred processing type (errorbar smaller than
        % symbol)
        plot( pos(sample_num)+posshift, r(sample_num,procA).(inputvar)(1), 'Marker',procMark ,'Color',colorlist{sample_num},'MarkerFaceColor',colorlist{sample_num},'markerSize',msA,'LineWidth',2*lw);
        hold on
        % EBSD values for individual scans (for "preferred processing type" only)
        plot( pos(sample_num)+posshift-smallshift, d(q(sample_num),procA).(inputvar), 'Marker',procMark,'Color',colorlist{sample_num},'markerSize',msA)
        hold on
        plot( pos(sample_num)+posshift+smallshift, d(q(sample_num)+1,procA).(inputvar), 'Marker',procMark,'Color',colorlist{sample_num},'markerSize',msA)
        hold on
        % EBSD data, non-preferred processing type 
        plot( pos(sample_num)+posshift, r(sample_num,procB).(inputvar)(1), 'Marker',procMarkB,'Color',colorlist{sample_num},'markerSize',msB,'LineWidth',2*lw)
        hold on
    end
    hold off
    xlim([xmin xmax])
    % ylim([0 4800])
    ylabel('d_T_m_t µm')
    set(gca,'xtick',[0 pos(1),pos(1)+posshift,pos(2),pos(2)+posshift],'xticklabel',names)
    ax=gca;
    ax.XTickLabelRotation = 90;
    ax.LineWidth = ax_lw;     
    ax.XAxis.FontSize = axx_fs;
    ax.YAxis.FontSize = axy_fs;

% Tmt Gbatch % % % % % % % % % % % % % % % % % % % % % % % % % % % 

    procA=1;
    procB=2;
    row=6;
    
    if procA==1
        procMark='x';
        procMarkB='o';
        msA=ms1;
        msB=ms2;
    else
        procMark='o';
        procMarkB='x';
        msA=ms2;
        msB=ms1;
    end
    
    inputvar='Gbatch_Tmt';
    
    subplot(2,2,4)
    for sample_num=1:2
        
        sel=TPM{sample_num};        
        
        % Pontesilli et al. data
        errorbar( pos(sample_num), sel(row,1), sel(row,2), sel(row,2) ,'Marker','s' ,'Color',colorlist{sample_num},'MarkerEdgeColor',colorlist{sample_num},'MarkerFaceColor',colorlist{sample_num},'markerSize',msS,'LineWidth',lw)
        hold on
        % EBSD data, preferred processing type (errorbar smaller than
        % symbol)
        plot( pos(sample_num)+posshift, 10*r(sample_num,procA).(inputvar)(1), 'Marker',procMark ,'Color',colorlist{sample_num},'MarkerFaceColor',colorlist{sample_num},'markerSize',msA,'LineWidth',2*lw);
        hold on
        % EBSD values for individual scans (for "preferred processing type" only)
        plot( pos(sample_num)+posshift-smallshift, 10*d(q(sample_num),procA).(inputvar), 'Marker',procMark,'Color',colorlist{sample_num},'markerSize',msA)
        hold on
        plot( pos(sample_num)+posshift+smallshift, 10*d(q(sample_num)+1,procA).(inputvar), 'Marker',procMark,'Color',colorlist{sample_num},'markerSize',msA)
        hold on
        % EBSD data, non-preferred processing type 
        plot( pos(sample_num)+posshift, 10*r(sample_num,procB).(inputvar)(1), 'Marker',procMarkB,'Color',colorlist{sample_num},'markerSize',msB,'LineWidth',2*lw)
        hold on
    end
    hold off
    xlim([xmin xmax])
    % ylim([0 4800])
    ylabel('G_b_a_t_c_h Tmt cms^-^1')
    set(gca,'xtick',[0 pos(1),pos(1)+posshift,pos(2),pos(2)+posshift],'xticklabel',names)
    ax=gca;
    ax.XTickLabelRotation = 90;
    ax.LineWidth = ax_lw;     
    ax.XAxis.FontSize = axx_fs;
    ax.YAxis.FontSize = axy_fs;

f = gcf;
f.Position = [leftfrac*w botfrac*h horzfrac*w vertfrac*0.631*h];

cd(savelocation)
exportgraphics(f,'fig7panel.png','Resolution',1000)
close
