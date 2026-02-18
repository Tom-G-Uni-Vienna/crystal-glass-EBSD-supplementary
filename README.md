# crystal-glass-EBSD-supplementary
Digital Supplementary Material to the paper "Optimizing Quantification of Microstructural Parameters of Crystal-Glass Assemblages by Electron Backscatter Diffraction", submitted for review to Journal of Petrology.

This repository contains Matlab scripts written using the MTEX toolbox version 5.11.2 (https://mtex-toolbox.github.io/) to assess the effect of different EBSD processing parameters on microstructural quantification of crystal-glass assembages from EBSD data. The repository also contains the raw EBSD data files used in the paper along with pre-generated Matlab variables (.m files) of the output of all processing scripts, allowing results to be immediately plotted with the plotting scripts.

## QUICK START GUIDE
To use the scripts provided you will need a copy of Matlab (unfortunately, Octave is not sufficient). MTEX recommends as recent a MATLAB version as possible, in any case newer than 2016b. The code was written & tested using version R2024b. Some scripts may require the “Statistics” MATLAB Toolbox. 

You then need to download and install the MTEX toolbox version 5.11.2 (instructions at https://mtex-toolbox.github.io/download). **Following installation, if you plan to run any of the scripts which re-import the raw data from .ang files, the standard MTEX version 5.11.2 should be modified by navigating to the folder _\mtex-5.11.2\interfaces\tools_ and replacing the file _loadHelper.m_ with the version provided in this repository.** This is required to prevent an error that can occur when importing .ang files, whereby not indexed points are not imported. The error is fixed in both MTEX versions 5.10 and 6, but present in 5.11.2.

If you are only interested in examining the processed results of the paper, initially by plotting them in the same way as the figures in the paper, you only need the scripts in folder ***MTEX scripts plots & tables***. These scripts use the MATLAB variables (outputs of the code) already saved in the folder MATLAB data. ***To speed up data import into MATLAB, make sure the entire folder structure is placed on an internal SSD hard disc!***

Only if you want to repeat the import of the .ang EBSD data files and the results of data processing do you need the scripts in folder ***MTEX scripts import and analyse***. A detailed explanation of the full sequence required to go from the raw data (folder Raw data) to the finished MATLAB data is found within that folder.

For any script, you will need to tell MATLAB where the relevant data is saved. This is accomplished by altering the very first variable in each script, supplementary_location, to be the location of the folder FINAL Supplementary on your system, following the format example used in the script.

## SCAN NAME GLOSSARY (INFO ALSO IN TABLE 1 OF MAIN PAPER)
*Original scan name on acquisition &rarr; Scan name in paper*\
ETK_5a_scan01 &rarr; Anhy_150\
ETK_5a_scan05 &rarr; Anhy_200\
ETK_5b_scan06 &rarr; Hy_200\
ETK_5b_scan07 &rarr; Hy_1000

## COMPLETE OUTLINE OF DIGITAL SUPPLEMENTARY MATERIAL
### Folder: Figures & tables
Empty folder, plotting scripts are set by default to output saved images to this folder, as is the script for generating supplementary table 1.

### Folder: MATLAB data
Contains all the necessary data for plotting all figures in the form of saved MATLAB variables (.mat files). These represent the pre-generated outputs of the scripts in the folder ***MTEX scripts import and analyse***, which can be directly used to generate all figure plots. You only need to run scripts in the folder ***MTEX scripts import and analyse*** if you want to repeat the process of generating the Matlab data.

### Folder: Measuring map distortion
Contains no MATLAB data. Here are saved the SE images of the 4 scanned areas which were used to measure the distortion of scans from their nominal sizes. Also included is a spreadsheet file with the measurements used. The explanation of which distances were measured is available in Supplementary Figure 1, which is also available in this folder for reference. 

### Folder: MTEX scripts import and analyse
This folder contains all scripts required to import the original .ang files into MATLAB/MTEX and distort them according to their true measured size and shape, as well as scripts designed to test the effect of different processing parameters and the scripts that generate the final “optimised” dataset. Data imported will be output into folder MATLAB data. A more detailed overview of what these scripts do and which scripts to use for what is available in the flowchart included in this folder. Each script itself is also extensively annotated. ***Please note that it is not necessary to carry out these scripts before running the plotting scripts (see quick start guide).***

### Folder: MTEX scripts plots & tables
This folder contains all scripts required to plot the results of data import and processing, including all figures and supplementary figures with the exception of supplementary figure 1 (which only contains the information on how to measure scan distortion).

### Folder: Raw data
This folder contains the EBSD scan data used in this study, exported from the OIM Analysis v. 8.0 software as .ang format (text) files. Files can be opened and viewed in Notepad and similar programs.

### File: loadHelper.m
**Following installation, if you plan to run any of the scripts which re-import the raw data from .ang files, the standard MTEX version 5.11.2 should be modified by navigating to the folder _\mtex-5.11.2\interfaces\tools_ and replacing the file _loadHelper.m_ with this version of the file.** This is required to prevent an error that can occur when importing .ang files, whereby not indexed points are not imported. The error is fixed in both MTEX versions 5.10 and 6.X, but present in 5.11.2.
