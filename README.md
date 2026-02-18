# crystal-glass-EBSD-supplementary
Digital Supplementary Material to the paper "Optimizing Quantification of Microstructural Parameters of Crystal-Glass Assemblages by Electron Backscatter Diffraction", submitted for review to Journal of Petrology.

This repository contains Matlab scripts written using the MTEX toolbox version 5.11.2 (https://mtex-toolbox.github.io/) to assess the effect of different EBSD processing parameters on microstructural quantification of crystal-glass assembages from EBSD data. The repository also contains the raw EBSD data files used as input and scripts that can be used to reproduce all plots used in the paper. 
> [!WARNING]
> ***The individual saved variables generated as outputs by these scripts are extremely large, some > 5 GB. If all scripts are run, the total size of the \MATLAB Data folder will be ca. 30 GB. Likewise, considerable computation time is required to generate each set of Matlab variables***

## HOW TO USE THESE SCRIPTS
To use the scripts provided you will need a copy of Matlab (unfortunately, Octave is not sufficient). MTEX recommends as recent a MATLAB version as possible, in any case newer than 2016b. The code was written & tested using version R2024b. Some scripts may require the “Statistics” MATLAB Toolbox. 

You then need to download and install the MTEX toolbox version 5.11.2 (instructions at https://mtex-toolbox.github.io/download). **Following installation, the standard MTEX version 5.11.2 should be modified by navigating to the folder _\mtex-5.11.2\interfaces\tools_ and replacing the file _loadHelper.m_ with the version provided in this repository.** This is required to prevent an error that can occur when importing .ang files, whereby not indexed points are not imported. The error is fixed in both MTEX versions 5.10 and 6, but present in 5.11.2.

> [!TIP]
> To speed up data import into MATLAB and subsequent saving of Matlab variables, make sure the entire folder structure is placed on an internal SSD hard disc!

Before running any script for the first time, you will need to tell MATLAB where the relevant data is saved. This is accomplished by altering the very first variable in each script, supplementary_location, to be the location of the folder ***\FINAL Supplementary*** on your system, following the format example used in the script.

Begin with the scripts in the folder ***\MTEX scripts import and analyse***. To import the raw data and obtain the final "optimised" microstructural quantification results, you only need to run script ***7_loop_import_regrid_and_clean_final.m*** followed by ***8_combine_scans.m***. Script 7 outputs three Matlab variable files, each ca. 1GB in size. The output of script 8 is a few hundred kb.

If you want to reproduce and plot the results of the full exploration of the effects of different EBSD processing parameters, you will also need to run the remaining scripts (numbered from 1 to 6) in the folder ***MTEX scripts import and analyse***. ***BEWARE, some of these scripts may take an hour or more to run and generate individual files of up to 5 GB in size as outputs. Make sure you have sufficient RAM and storage space on your computer before running these scripts***. However, each script must only be run once and its output saved to disk.

A flowchart indicating which of scripts 1 to 6 must be run in which order to generate which final results is included in this repository as file ***Import_and_analyse_scripts_flowchart*** (in both .pdf and .odp format)

In addition, if you want to generate the images used in figures 1, 8, 9, S5 and S6, you must run the scripts contained in the subfolder ***\Import undistorted*** in order from A to C. The only difference with these scripts is that they do not distort the EBSD maps to their true size and shape. Therefore, the outputs retain hexagonal pixels which make them better suited for illustrative figures.

In folder ***\MTEX scripts plots & tables*** you can find all scripts required to plot the results of data import and processing, including all figures and supplementary figures with the exception of supplementary figure 1 (which only contains the information on how to measure scan distortion), and supplementary table 1. Please check the beginning of each plotting script to see which Matlab variables must be generated before plotting.

## COMPLETE OUTLINE OF DIGITAL SUPPLEMENTARY MATERIAL
### Folder: Figures & tables
Empty folder, plotting scripts are set by default to output saved images to this folder, as is the script for generating supplementary table 1.

### Folder: MATLAB data
Empty folder, here the scripts in the folder ***MTEX scripts import and analyse*** will save the necessary data for plotting all figures in the form of MATLAB variables (.mat files).

### Folder: MTEX scripts import and analyse
This folder contains all scripts required to import the original .ang files into MATLAB/MTEX and distort them according to their true measured size and shape, as well as scripts designed to test the effect of different processing parameters and the scripts that generate the final “optimised” dataset. Data imported will be output into folder ***\MATLAB data***. A more detailed overview of what these scripts do and which scripts to use for what is available in the flowchart included in this repository. Each script itself is also extensively annotated.

### Folder: MTEX scripts plots & tables
This folder contains all scripts required to plot the results of data import and processing, including all figures and supplementary figures with the exception of supplementary figure 1 (which only contains the information on how to measure scan distortion).

### Folder: Raw data
This folder contains the EBSD scan data used in this study, exported from the OIM Analysis v. 8.0 software as .ang format (text) files and compressed as .7z files. Following decompression, files can be opened and viewed in Notepad and similar programs. Note that for other scripts to find these files, they must be place directly in this folders, not in subfolders.

### File: loadHelper.m
**Following installation, the standard MTEX version 5.11.2 should be modified by navigating to the folder _\mtex-5.11.2\interfaces\tools_ and replacing the file _loadHelper.m_ with this version of the file.** This is required to prevent an error that can occur when importing .ang files, whereby not indexed points are not imported. The error is fixed in both MTEX versions 5.10 and 6.X, but present in 5.11.2.

## SCAN NAME GLOSSARY (INFO ALSO IN TABLE 1 OF MAIN PAPER)
*Original scan name on acquisition &rarr; Scan name in paper*\
ETK_5a_scan01 &rarr; Anhy_150\
ETK_5a_scan05 &rarr; Anhy_200\
ETK_5b_scan06 &rarr; Hy_200\
ETK_5b_scan07 &rarr; Hy_1000
