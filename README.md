# Seahorse-Analysis
Analysis of metabolic profiles and calculation of parameters from output of multiple Seahorse experiments.

SHorE: Seahorse Explorer                                                                                                                
Author: Elizabeth A. Proctor                                                 
proctor.tools@gmail.com 
GitHub: elizabethproctor 

This software is protected under the Gnu Public License, version 3 (GPLv3) 
**********
Please cite as:
D. Nicholas, E. A. Proctor, F. Raval, B. C. Ip, C. Habib, E. Ritou, T. N. Grammatopoulos, C. M. Apovian, D. A. Lauffenburger, B. N. Nikolajczyk. "Advances in the Quantification of Mitochondrial Function in Primary Human Immune Cells through Extracellular Flux Analysis," submitted (2016)
**********

This program combines results from multiple experiments performed by Seahorse XF Cell Mitochondrial Stress Test, and calculates unified OCR:ECAR ratios and metabolic parameters for each sample/condition:
	(1) Basal Respiration (BR), OCR:ECAR for Basal Respiration
	(2) ATP Production (ATP)
	(3) Maximal Respiration (MR), OCR:ECAR for Maximal Respiration
	(4) Spare Respiratory Capacity (SRC)
	(5) Proton Leak (PL)
	(6) Non-Mitochondrial Respiration (NMiR)


To use this program, you will need:
(1) For each Seahorse experiment, the "Rate" tab of Seahorse output exported in CSV format. Add a column with corresponding cell numbers plated (or other normalization information you wish to use) for each well.
(2) A single text file listing the filenames of each CSV file from 1 above.


Instructions for using command line:
MAC/LINUX/UNIX
	1)	Open the Terminal program. In recent Mac OSX versions this application is under the “Utilities” sub-folder in Applications.
	2)	If you are in the folder with your data (denoted by the file path that appears next to your prompt), skip this step. If not, type “cd “ followed by the path to where your data is. For example:
		cd ~/projects/seahorse
	3)	Follow instructions for “Command line usage,” below. Example command line, for when the program is located in the same directory as all data files:
		./SHorE.py
	followed by the arguments specific to your data, as described below

WINDOWS
	1)	Open Windows “command prompt”: Start > Programs > Accessories > Command Prompt. You should see something like “C:\>”, or similar.
	2)	Type:
		py
	and hit enter. If you get an error, you will need to download and install Python 2, which is free: https://www.python.org/downloads/. If no error occurs, you have Python already installed. Type ctrl+D to exit.
	3)	If you are in the folder with your data (denoted by the file path that appears next to your prompt), skip this step. If not, type “cd “ followed by the path to where your data is. For example:
		cd \Users\username\projects\seahorse
	4)	Follow instructions for “Command line usage,” below. Example command line, for when the program is located in the same directory as all data files:
		py SHorE.py
	

Command line usage: SHorE.py list outputName(optional) nMeasure(optional)
 	"list" is the text file described in “To use this program (2),” above.
	"outputName" is a root for the name of the output files that will be generated. This argument is optional, and is not required for the program to run. Default is "log.out," with NO OVERWRITE PROTECTION.
	“nMeasure” is the total number of measurements performed. This argument is optional, and is not required for the program to run. Default is 14, laid out as below. 12 is also implemented, but needs to be specified.

NOTE: CONDITION NAMES THAT GO WITH CELL COUNTS MUST BE WRITTEN EXACTLY AS IN THE SEAHORSE OUTPUT FILE

NOTE: THERE IS NO OVERWRITE PROTECTION FOR OUTPUT FILES. If you do not list a different output name and run the program twice in the same directory, your first output files will be overwritten.


Program notes
This program assumes the following experimental layout (14-measurement format):
	Measurement time points, separated by 5 minutes:
		1-5: Baseline
		---> Inject oligomycin (ATP coupler)
		6-8
		---> Inject FCCP (ETC accelerator)
		9-11
		---> Inject antimycin/rotenone (mitochondrial inhibitor)
		12-14
12-measurement format is similar, with two fewer baseline measurements.

Handling of measurement error:
	⁃	Errors are calculated as standard errors of the mean, with commonly used error propogation for calculated parameters. 
	⁃	Final oxygen consumption curves are determined as medians of biological and technical replicates, in order to discard outliers and account for potential skewed distributions. In a perfectly normal distribution, median = mean.
	⁃	For metabolic parameters in which a median of similar measurements is used (e.g. NMiR), the error for the selected measurement is designated as the error for the parameter.

For comparison purposes, each replicate across all experiments is aligned to the non-metabolic respiration of the first instance of that condition. Measurements are not otherwise scaled.
For each condition, each measurement is represented by the median of replicates, across all experiments.
Metabolic parameters are calculated as:
	•	Non-Mitochondrial Respiration (NMiR) is calculated as the median of measurements 12-14, following injection of mitochondrial inhibitor. The median is used in order to discount potential non-sense measurements due to the low level of respiration.
	•	Basal Respiration (BR) is calculated as average of timepoints 3-5, minus Non-Mitochondrial Respiration (NMiR). For 12 measurements, BR is the average of timepoints 2 and 3.
	•	ATP Production (ATP) is calculated as the difference between Basal Respiration and the median of measurements 6-8, which are taken following injection of ATP coupler. The median is used in order to discount potential non-sense measurements due to the low level of respiration.
	•	Maximal Respiration (MR) is calculated as the first measurement after injection of ETC accelerator (timepoint 9), minus Non-Mitochondrial Respiration (NMiR).
	•	Spare Respiratory Capacity (SRC) is calculated as the difference between Maximal Respiration (MR) and the average of measurements 3-5, which represent the baseline (or Basal Respiration + Non-Mitochondrial Respiration).
Wells with OCR or ECAR less than or equal to zero (at any measurement) will be excluded from calculations.
