#!/usr/bin/python

####################################################################################
##   SHorE: Seahorse Explorer                                                     ##
##                                                                                ##
##           Author: Elizabeth A. Proctor                                         ##
##           proctor.tools@gmail.com                                              ##
##           GitHub: elizabethproctor                                             ##
##                                                                                ##
##   Please cite as:                                                              ##
##           D. Nicholas, E. A. Proctor, F. M. Raval, B. C. Ip, C. Habib,         ##
##           E. Ritou, T. N. Grammatopoulos, D. Steenkamp, H. Dooms,              ##
##           C. M. Apovian, D. A. Lauffenburger, B. S. Nikolajczyk. "Advances in  ##
##           the quantification of mitochondrial function in primary human immune ##
##           cells through extracellular flux analysis," PLoS ONE 12(2):e0170975  ##
##           (2017)                                                               ##
##                                                                                ##
##   This software is protected under the Gnu Public License, version 3 (GPLv3)   ##
####################################################################################


# This program combines results from multiple experiments performed by Seahorse XF Cell Mitochondrial Stress Test, and calculates unified OCR:ECAR ratios and metabolic parameters for each sample/condition:
#			(1) Basal Respiration (BR), OCR:ECAR for Basal Respiration
#   	(2) ATP Production (ATP)
#   	(3) Maximal Respiration (MR), OCR:ECAR for Maximal Respiration
#   	(4) Spare Respiratory Capacity (SRC)
#   	(5) Proton Leak (PL)
#   	(6) Non-Mitochondrial Respiration (NMiR)
#
# To use this program, you will need:
#			(1) For each Seahorse experiment, the "Rate" tab of Seahorse output exported in CSV format. Add a column with corresponding cell numbers plated (or other normalization information you wish to use) for each well.
#			(2) A single text file listing the filenames of each CSV file from #1 above.
#
# Command line usage: SHorE.py list outputName(optional) nMeasure(optional)
# 	"list" is the text file described in (2), above.
#		"outputName" is the name of the output file that will be generated. This argument is optional, and is not required for the program to run. Default is "log.out," with NO OVERWRITE PROTECTION.
#		nMeasure is the total number of measurements performed. This argument is optional, and is not required for the program to run. Default is 14, laid out as below. 12 is also implemented, but needs to be specified.
#		NOTE: CONDITION NAMES THAT GO WITH CELL COUNTS MUST BE WRITTEN EXACTLY AS IN THE SEAHORSE OUTPUT FILE
#
# This program assumes the following experimental layout (14-measurement format):
# 	Measurement time points, separated by 5 minutes:
#			1-5: Baseline
#			---> Inject oligomycin (ATP coupler)
#			6-8
#			---> Inject FCCP (ETC accelerator)
#			9-11
#			---> Inject antimycin/rotenone (mitochondrial inhibitor)
#			12-14
#		12-measurement format has two fewer baseline measurements.
#
# Handling of measurement error:
# 	- Errors are calculated as standard errors of the mean, with commonly used error propogation for calculated parameters. 
# 	- Final oxygen consumption curves are determined as medians of biological and technical replicates, in order to discard outliers and account for potential skewed distributions. In a perfectly normal distribution, median = mean.
# 	- For metabolic parameters in which a median of similar measurements is used (e.g. NMiR), the error for the selected measurement is designated as the error for the parameter.
#
###################
#	PROGRAM NOTES 
###################
# For comparison purposes, each replicate across all experiments is aligned to the non-metabolic respiration of the first instance of that condition. Measurements are not otherwise scaled.
# For each condition, each measurement is represented by the median of replicates, across all experiments.
# Metabolic parameters are calculated as:
#		Non-Mitochondrial Respiration (NMiR) is calculated as the median of measurements 12-14, following injection of mitochondrial inhibitor. 
#				The median is used in order to discount potential non-sense measurements due to the low level of respiration.
#		Basal Respiration (BR) is calculated as average of timepoints 3-5, minus Non-Mitochondrial Respiration (NMiR). For 12 measurements, BR is the average of timepoints 2 and 3.
# 	ATP Production (ATP) is calculated as the difference between Basal Respiration and the median of measurements 6-8, which are taken following injection of ATP coupler.
#				The median is used in order to discount potential non-sense measurements due to the low level of respiration.
# 	Maximal Respiration (MR) is calculated as the first measurement after injection of ETC accelerator (timepoint 9), minus Non-Mitochondrial Respiration (NMiR).
#		Spare Respiratory Capacity (SRC) is calculated as the difference between Maximal Respiration (MR) and the average of measurements 3-5, which represent the baseline (or Basal Respiration + Non-Mitochondrial Respiration).
# Wells with OCR or ECAR less than or equal to zero (at any measurement) will be excluded from calculations.

from sys import argv, exit
import csv
import numpy
from decimal import Decimal

def median(a):
	sorta = sorted(a)
	length = len(a)
	if length % 2:
		ans = sorta[(length-1)/2]
	else:
		ans = (sorta[(length/2)-1]+sorta[length/2])/2
	return ans

def mean(a):
	ans = sum(a)/float(len(a))
	return ans

# Create list of conditions
conditions = []
try:
	f = open(argv[1],'r').readlines()  # Open the input file listing filenames of CSV files from individual Seahorse experiments
except:
	exit("ERROR: No input file (list of Seahorse CSV filenames) given.")
if len(argv) > 2:
	logname = str(argv[2])
else:
	logname = "log.out"
wellLog = "wells_"+logname
groupLog = "groups_"+logname
gw = open(wellLog,'w')	# Create well-by-well output file
gg = open(groupLog,'w')  # Create group output file
for line in f:
	# Conditions are read for EACH file and inserted into list of conditions before the rest of the program begins, so that there is only one instance of each condition present in all experiments combined.
	split = line.split()
	try:
		fname = split[0]
	except:
		exit("ERROR: Faulty filename given in input file. Check for and remove empty lines.")
	with open(fname, 'rU') as g:
		reader = csv.reader(g)
		for row in reader:
			if str(row[2]) not in conditions and str(row[2])!="Group":
				conditions.append(str(row[2]))
print "Total number of conditions:", len(conditions)

# Define measurement regime
nMeas = 14
if len(argv) > 3:
	nMeas = int(argv[3])

# Create matrices to hold data read-in  
reference = [[] for j in range(len(conditions))] # initialize matrix for reference (N conditions x M measurements), for baselining between files.
reference2 = [[] for j in range(len(conditions))] # initialize matrix for ECAR reference
working = [[[] for i in range(nMeas)] for j in range(len(conditions))] # initialize matrix for working record of final OCRs for each condition/measurement for each file after being baselined to reference (N conditions x M measurements x P replicates) 
working2 = [[[] for i in range(nMeas)] for j in range(len(conditions))]
nonnegConsumpRatios = [[[] for i in range(nMeas)] for j in range(len(conditions))]
gw = open(wellLog,'a')
gg = open(groupLog,'a')

# Start reading experiment files
for line in f:	# For each line in input file naming experiment files, i.e. for each experiment
	split = line.split()
	fname = split[0].rstrip()
	normaveOCR = [[] for j in range(len(conditions))] # initialize matrix for N conditions, for which vector of replicates for each measurement will be appended after normalizing for cell count.
	normaveECAR = [[] for j in range(len(conditions))] # initialize matrix for N conditions, for which vector of replicates for each measurement will be appended after normalizing for cell count.
	allConsumpRatios = [[[] for i in range(nMeas)] for j in range(len(conditions))] # initialize matrix for N conditions x M measurements
	# Initialize temporary 3D matrix of measurement x condition x replicate. Note non-constant number of replicates.
	temp = [[[] for i in range(nMeas)] for j in range(len(conditions))]	# Create matrix "temp", which will house the OCR normalized by cell count for each individual experiment. Note that this matrix is redefined for each iteration of this loop (each individual experiment), so even though elements are added by "append", the matrix is cleared at the end of the experiment read.
	temp2 = [[[] for i in range(nMeas)] for j in range(len(conditions))] # Create same matrix, but for ECAR, for use later calculating OCR:ECAR ratio
	with open(fname, 'rU') as h: # read CSV file
		reader = csv.reader(h)
		for row in reader:
			if str(row[0]) != "Measurement": # skip the first line, which is just headers and not experimental data
				try:
					measurement = int(row[0]) # measurement number
				except:
					exit("ERROR: CSV file may be empty or not contain all measurements. Please check that columns Measurement, Well, Group, Time, OCR, ECAR, PPR, Cell Count are present and in order.")
				try:
					condition = str(row[2])
				except:
					exit("ERROR: CSV file does not contain all measurements. Please check that columns Measurement, Well, Group, Time, OCR, ECAR, PPR, Cell Count are present and in order.")
				try:
					ocr = float(row[4])
				except:
					exit("ERROR: CSV file does not contain all measurements. Please check that columns Measurement, Well, Group, Time, OCR, ECAR, PPR, Cell Count are present and in order.")
				try:
					ecar = float(row[5])
				except:
					exit("ERROR: CSV file does not contain all measurements. Please check that columns Measurement, Well, Group, Time, OCR, ECAR, PPR, Cell Count are present and in order.")
				try:
					if row[7] == "":
						cellnum = 0
					else:
						cellnum = float(row[7])	# Grab cell counts for baseline from last column of CSV
				except:
					exit("ERROR: Please include cell counts for each well in the eighth column of each CSV file. Cell numbers are used for normalization of oxygen consumption measurements.")
				if (cellnum != 0) and (ocr >= 0) and (ecar > 0):
					temp[conditions.index(condition)][measurement-1].append(ocr/cellnum)	# At once normalize OCR by cell count and insert into temporary matrix.
					temp2[conditions.index(condition)][measurement-1].append(ecar/cellnum)	# Same operations for creation of ECAR matrix
					allConsumpRatios[conditions.index(condition)][measurement-1].append(ocr/ecar)
				elif cellnum == 0:
					temp[conditions.index(condition)][measurement-1].append(0)	# If corresponding cell count is zero, then OCR should automatically also be zero.
					temp2[conditions.index(condition)][measurement-1].append(0)
					allConsumpRatios[conditions.index(condition)][measurement-1].append(0)
				else:
					temp[conditions.index(condition)][measurement-1].append(-1)
					temp2[conditions.index(condition)][measurement-1].append(-1)
					allConsumpRatios[conditions.index(condition)][measurement-1].append(-1)
	for i in range(len(conditions)):
		if (conditions[i] != "Background") and (conditions[i] != "Unassigned"):
			badOCR = []
			for j in range(nMeas):
				for k in range(len(temp[i][j])):
					if temp[i][j][k] == -1:
						badOCR.append(k)
					else:
						pass
			if badOCR != []:
				badOCR_noDups = []
				for n in badOCR:
					if n not in badOCR_noDups:
						badOCR_noDups.append(n)
				badOCR_noDups.sort()
				count = 0
				for m in badOCR_noDups:
					for j in range(nMeas):
						del temp[i][j][m-count]
						del temp2[i][j][m-count]
						del allConsumpRatios[i][j][m-count] 
					count += 1
	# Filter out empty lists into new OCR array, calculate oxygen consumption ratio
	for i in range(len(conditions)):
		for j in range(nMeas):
			if allConsumpRatios[i][j] != []:
				for k in allConsumpRatios[i][j]:
					nonnegConsumpRatios[i][j].append(k)
			else:
				pass
			if temp[i][j] != []:
				normaveOCR[i].append(temp[i][j])
				normaveECAR[i].append(temp2[i][j])
			else:
				pass
	# Baseline by Non-Mitochondrial Respiration of first biological replicate for each group/condition
	for i in range(len(conditions)):
		if (len(normaveOCR[i]) < nMeas) and (normaveOCR[i] != []):
			print "SKIPPING CONDITION", conditions[i], "IN EXPERIMENT", line.rstrip(), ": INCOMPLETE MEASUREMENTS"
			continue
		if (reference[i] == []) and (temp[i][0] != []) and (normaveOCR[i] != []):	# If this is the first experiment in the loop, make the first replicate for this condition the reference for baselining.
			for j in range(nMeas):	# Append the values for this experiment to the "working" matrix, which will store normalized, baselined values for all experiments, which will be combined later by median.
				reference[i].append(normaveOCR[i][j][0])
				reference2[i].append(normaveECAR[i][j][0])
			if nMeas == 14:
				refNMiR = median([reference[i][11], reference[i][12], reference[i][13]])
				for k in range(len(normaveOCR[i][0])):
					instNMiR = median([normaveOCR[i][11][k], normaveOCR[i][12][k], normaveOCR[i][13][k]])
					diff = instNMiR - refNMiR
					for j in range(nMeas):
						working[i][j].append(normaveOCR[i][j][k] - diff)
						working2[i][j].append(normaveECAR[i][j][k])
			elif nMeas == 12:
				refNMiR = median([reference[i][9], reference[i][10], reference[i][11]])
				for k in range(len(normaveOCR[i][0])):
					instNMiR = median([normaveOCR[i][9][k], normaveOCR[i][10][k], normaveOCR[i][11][k]])
					diff = instNMiR - refNMiR
					for j in range(nMeas):
						working[i][j].append(normaveOCR[i][j][k] - diff)
						working2[i][j].append(normaveECAR[i][j][k])
			else:
				print "ERROR: incorrect number of measurements specified. Only 14 and 12 are currently implemented."
				raise SystemExit
		else: # if a previous experiment file already had data for this condition
			# Set non-mitochondrial respiration to equal that of the reference. Record difference, and add/subtract from all points for that condition.
			if nMeas == 14:
				if normaveOCR[i] != []:
					refNMiR = median([reference[i][11], reference[i][12], reference[i][13]])
					for k in range(len(normaveOCR[i][0])):
						instNMiR = median([normaveOCR[i][11][k], normaveOCR[i][12][k], normaveOCR[i][13][k]])
						diff = instNMiR - refNMiR
						diff2 = normaveECAR[i][j][k] - reference2[i][j]
						# Insert baselined points into "working" matrix for that condition.
						for j in range(nMeas):
							working[i][j].append(normaveOCR[i][j][k]-diff)
							working2[i][j].append(normaveECAR[i][j][k]-diff2)
				else:
					pass
			elif nMeas == 12:
				if normaveOCR[i] != []:
					refNMiR = median([reference[i][9], reference[i][10], reference[i][11]])
					for k in range(len(normaveOCR[i][0])):
						instNMiR = median([normaveOCR[i][9][k], normaveOCR[i][10][k], normaveOCR[i][11][k]])
						diff = instNMiR - refNMiR
						diff2 = normaveECAR[i][j][k] - reference2[i][j][k]
						for j in range(nMeas):
							working[i][j].append(normaveOCR[i][j][k]-diff)
							working2[i][j].append(normaveECAR[i][j][k]-diff2)
				else:
					pass
			else:
				print "ERROR: incorrect number of measurements. Only 14 and 12 are currently implemented."
				raise SystemExit

finalConsumpRatios = [[0 for i in range(nMeas)] for j in range(len(conditions))]
consumpRatioErr = [[0 for i in range(nMeas)] for j in range(len(conditions))]
for i in range(len(conditions)):
	for j in range(nMeas):
		if nonnegConsumpRatios[i][j] == []:
			pass
		else:
			finalConsumpRatios[i][j] = median(nonnegConsumpRatios[i][j])
			consumpRatioErr[i][j] = numpy.std(nonnegConsumpRatios[i][j])/numpy.sqrt(len(nonnegConsumpRatios[i][j]))

# Take median of each measurement across "working" to get final of that measurement. Error is calculated as SEM.
finalM = [[0 for i in range(nMeas)] for j in range(len(conditions))] # initialize matrix for median of measurements for each condition across files
OCRerr = [[0 for i in range(nMeas)] for j in range(len(conditions))]
finalMe = [[0 for i in range(nMeas)] for j in range(len(conditions))]
ECARerr = [[0 for i in range(nMeas)] for j in range(len(conditions))]
gw.write("##########"+'\n'+"WELL-BY-WELL MEASUREMENT CURVES"+'\n'+"##########"+'\n')
for i in range(len(conditions)):
	if (conditions[i] != "Background") and (conditions[i] != "Unassigned"):
		gw.write('\n'+"Condition: "+conditions[i]+'\n')
		if working[i][0] == []:
			gw.write('\t'+"No complete curves for this condition: all replicates contain at least one negative OCR measurement. Please check input data."+'\n')
		for k in range(len(working[i][0])):
			gw.write('\t'+"Replicate "+str(k+1)+'\n'+'\t'+"## OCR ##"+'\n')
			for j in range(nMeas):
				gw.write('\t'+'\t'+"Measurement "+str(j+1)+": "+str(working[i][j][k])+'\n')
				finalM[i][j] = median(working[i][j])
				OCRerr[i][j] = numpy.std(working[i][j])/numpy.sqrt(len(working[i][j]))
			gw.write('\t'+"## ECAR ##"+'\n')
			for j in range(nMeas):
				gw.write('\t'+'\t'+"Measurement "+str(j+1)+": "+str(working2[i][j][k])+'\n')
				finalMe[i][j] = median(working2[i][j])
				ECARerr[i][j] = numpy.std(working[i][j])/numpy.sqrt(len(working2[i][j]))
gw.write('\n'+'\n'+"##########"+'\n'+"WELL-BY-WELL RESPIRATION CALCULATIONS"+'\n'+"##########"+'\n')

# Calculate well-by-well metabolic parameters
for i in range(len(conditions)):
	if (conditions[i] != "Background") and (conditions[i] != "Unassigned"):
		gw.write('\n'+"Condition: "+conditions[i]+'\n')
		if working[i][0] == []:
			gw.write('\t'+'\t'+"No metabolic parameters could be calculated because all replicates contain at least one negative OCR measurement. Please check input data.")
		for k in range(len(working[i][0])):
			gw.write('\t'+"Replicate "+str(k+1)+'\n')
			if nMeas == 14:
				# Calculate Non-Mitorchondrial Respiration (well-by-well)
				NMiRi = median([working[i][11][k], working[i][12][k], working[i][13][k]])
				if NMiRi<0:
					NMiRi = 0
				gw.write('\t'+'\t'+"Non-Mitochondrial Respiraion = "+str(NMiRi)+'\n')
				# Calculate Basal Respiration (well-by-well)
				BRi = mean([working[i][2][k], working[i][3][k], working[i][4][k]]) - NMiRi
				gw.write('\t'+'\t'+"Basal Respiration = "+str(BRi)+'\n')
				gw.write('\t'+'\t'+'\t'+"OCR:ECAR = "+str(mean([nonnegConsumpRatios[i][2][k], nonnegConsumpRatios[i][3][k], nonnegConsumpRatios[i][4][k]]))+'\n')
				# Calculate Proton Leak (well-by-well)
				PLi = median([working[i][5][k], working[i][6][k], working[i][7][k]])
				if PLi<0:
					PLi = 0
				gw.write('\t'+'\t'+"Proton Leak = "+str(PLi)+'\n')
				# Calculate ATP Production (well-by-well)
				ATPi = BRi - PLi
				gw.write('\t'+'\t'+"ATP Production = "+str(ATPi)+'\n')
				# Calculate Maximal Respiration
				MRi = max(working[i][8][k], working[i][9][k], working[i][10][k]) - NMiRi
				gw.write('\t'+'\t'+"Maximal Respiration = "+str(MRi)+'\n')
				junk = max(working[i][8][k], working[i][9][k], working[i][10][k])
				for m in [8, 9, 10]:
					if junk == working[i][m][k]:
						gw.write('\t'+'\t'+'\t'+"OCR:ECAR = "+str(nonnegConsumpRatios[i][m][k])+'\n')
				# Sanity check
				if BRi > MRi:
					print "ERROR: Maximal Respiration is lower than Basal Respiration!", conditions[i], "Replicate", str(k+1)
				# Calculate Spare Respiratory Capacity
				SRCi = MRi - BRi
				gw.write('\t'+'\t'+"Spare Respiratory Capacity = "+str(SRCi)+'\n')
			elif nMeas == 12:
				# Calculate Non-Mitorchondrial Respiration (well-by-well)
				NMiRi = median([working[i][9][k], working[i][10][k], working[i][11][k]])
				if NMiRi<0:
					NMiRi = 0
				gw.write('\t'+'\t'+"Non-Mitochondrial Respiraion = "+str(NMiRi)+'\n')
				# Calculate Basal Respiration (well-by-well)
				BRi = mean([working[i][1][k], working[i][2][k]]) - NMiRi
				gw.write('\t'+'\t'+"Basal Respiration = "+str(BRi)+'\n')
				gw.write('\t'+'\t'+'\t'+"OCR:ECAR = "+str(mean([nonnegConsumpRatios[i][1][k], nonnegConsumpRatios[i][2][k]]))+'\n')
				# Calculate Proton Leak (well-by-well)
				PLi = median([working[i][3][k], working[i][4][k], working[i][5][k]])
				if PLi<0:
					PLi = 0
				gw.write('\t'+'\t'+"Proton Leak = "+str(PLi)+'\n')
				# Calculate ATP Production (well-by-well)
				ATPi = BRi - PLi
				gw.write('\t'+'\t'+"ATP Production = "+str(ATPi)+'\n')
				# Calculate Maximal Respiration
				MRi = max(working[i][6][k], working[i][7][k], working[i][8][k]) - NMiRi
				gw.write('\t'+'\t'+"Maximal Respiration = "+str(MRi)+'\n')
				junk = max(working[i][6][k], working[i][7][k], working[i][8][k])
				for m in [6, 7, 8]:
					if junk == working[i][m][k]:
						gw.write('\t'+'\t'+'\t'+"OCR:ECAR = "+str(nonnegConsumpRatios[i][m][k])+'\n')
				# Sanity check
				if BRi > MRi:
					print "ERROR: Maximal Respiration is lower than Basal Respiration!", conditions[i], "Replicate", str(k+1)
				# Calculate Spare Respiratory Capacity
				SRCi = MRi - BRi
				gw.write('\t'+'\t'+"Spare Respiratory Capacity = "+str(SRCi)+'\n')


# Record final group oxygen consumption ratios
gg.write("##########"+'\n'+"GROUP MEASUREMENT CURVES"+'\n'+"##########"+'\n')
for i in range(len(conditions)):
	if (conditions[i] != "Background") and (conditions[i] != "Unassigned"):
		gg.write('\n'+conditions[i]+": "+'\n'+"## OCR #"+'\n')
		for j in range(nMeas):
			gg.write('\t'+"Measurement "+str(j+1)+": "+str(finalM[i][j])+" +- "+str(OCRerr[i][j])+'\n')
		gg.write("## ECAR ##"+'\n')
		for j in range(nMeas):
			gg.write('\t'+"Measurement "+str(j+1)+": "+str(finalMe[i][j])+" +- "+str(ECARerr[i][j])+'\n')

gg.write('\n'+'\n'+"##########"+'\n'+"GROUP RESPIRATION CALCULATIONS"+'\n'+"##########"+'\n')

# Final respiratory calculations
for i in range(len(conditions)):
	gg.write('\n'+conditions[i]+'\n')
	if nMeas ==14:
		# Calculate Non-Mitochondrial Respiration
		NMiR = median([finalM[i][11], finalM[i][12], finalM[i][13]])
		NMiRerr = OCRerr[i][finalM[i].index(NMiR)]
		if NMiR<0:
			NMiR = 0
		gg.write('\t'+"Non-Mitochondrial Respiration = "+str(NMiR)+" +- "+str(NMiRerr)+'\n')
		# Calculate Basal Respiration
		BR = mean([finalM[i][2], finalM[i][3], finalM[i][4]]) - NMiR
		BRerr = numpy.sqrt((numpy.sqrt(OCRerr[i][2]**2+OCRerr[i][3]**2+OCRerr[i][4]**2)/3)**2 + NMiRerr**2)
		gg.write('\t'+"Basal Respiration = "+str(BR)+" +- "+str(BRerr)+'\n')
		gg.write('\t'+'\t'+"OCR:ECAR = "+str(mean([finalConsumpRatios[i][2], finalConsumpRatios[i][3], finalConsumpRatios[i][4]]))+" +- "+str(mean([consumpRatioErr[i][2], consumpRatioErr[i][3], consumpRatioErr[i][4]]))+'\n')
		# Calculate Proton Leak
		PL = median([finalM[i][5], finalM[i][6], finalM[i][7]]) - NMiR
		junk = median([finalM[i][5], finalM[i][6], finalM[i][7]])
		PLerr = numpy.sqrt(OCRerr[i][finalM[i].index(junk)]**2 + NMiRerr**2)
		if PL<0:
			PL = 0
		gg.write('\t'+"Proton Leak = "+str(PL)+" +- "+str(PLerr)+'\n')
		# Calculate ATP Production
		ATP = BR - PL
		ATPerr = numpy.sqrt(BRerr**2 + PLerr**2)
		gg.write('\t'+"ATP Production = "+str(ATP)+" +- "+str(ATPerr)+'\n')
		# Calculate Maximal Respiration
		MR = max(finalM[i][8], finalM[i][9], finalM[i][10]) - NMiR
		junk = max(finalM[i][8], finalM[i][9], finalM[i][10])
		MRerr = numpy.sqrt(OCRerr[i][finalM[i].index(junk)]**2 + NMiRerr**2)
		gg.write('\t'+"Maximal Respiration = "+str(MR)+" +- "+str(MRerr)+'\n')
		gg.write('\t'+'\t'+"OCR:ECAR = "+str(finalConsumpRatios[i][finalM[i].index(junk)])+" +- "+str(consumpRatioErr[i][finalM[i].index(junk)])+'\n')
		# Sanity check
		if BR > MR:
			print "ERROR: Maximal Respiration is lower than Basal Respiration!", conditions[i]
		# Calculate Spare Respiratory Capacity
		SRC = MR - BR
		SRCerr = numpy.sqrt(MRerr**2 + BRerr**2)
		gg.write('\t'+"Spare Respiratory Capacity = "+str(SRC)+" +- "+str(SRCerr)+'\n')
	elif nMeas == 12:
		# Calculate Non-Mitochondrial Respiration
		NMiR = median([finalM[i][9], finalM[i][10], finalM[i][11]])
		NMiRerr = OCRerr[i][finalM[i].index(NMiR)]
		if NMiR<0:
		  NMiR = 0
		gg.write('\t'+"Non-Mitochondrial Respiration = "+str(NMiR)+" +- "+str(NMiRerr)+'\n')
		# Calculate Basal Respiration
		BR = mean([finalM[i][1], finalM[i][2]]) - NMiR
		BRerr = numpy.sqrt((numpy.sqrt(OCRerr[i][1]**2+OCRerr[i][2]**2)/2)**2 + NMiRerr**2)
		gg.write('\t'+"Basal Respiration = "+str(BR)+" +- "+str(BRerr)+'\n')
		gg.write('\t'+'\t'+"OCR:ECAR = "+str(mean([finalConsumpRatios[i][1], finalConsumpRatios[i][2]]))+" +- "+str(mean([consumpRatioErr[i][1], consumpRatioErr[i][2]]))+'\n')
		# Calculate Proton Leak
		PL = median([finalM[i][3], finalM[i][4], finalM[i][5]]) - NMiR
		junk = median([finalM[i][3], finalM[i][4], finalM[i][5]])
		PLerr = numpy.sqrt(OCRerr[i][finalM[i].index(junk)]**2 + NMiRerr**2)
		if PL<0:
		  PL = 0
		gg.write('\t'+"Proton Leak = "+str(PL)+" +- "+str(PLerr)+'\n')
		# Calculate ATP Production
		ATP = BR - PL
		ATPerr = numpy.sqrt(BRerr**2 + PLerr**2)
		gg.write('\t'+"ATP Production = "+str(ATP)+" +- "+str(ATPerr)+'\n')
		# Calculate Maximal Respiration
		MR = max(finalM[i][6], finalM[i][7], finalM[i][8]) - NMiR
		junk = max(finalM[i][6], finalM[i][7], finalM[i][8])
		MRerr = numpy.sqrt(OCRerr[i][finalM[i].index(junk)]**2 + NMiRerr**2)
		gg.write('\t'+"Maximal Respiration = "+str(MR)+" +- "+str(MRerr)+'\n')
		gg.write('\t'+'\t'+"OCR:ECAR = "+str(finalConsumpRatios[i][finalM[i].index(junk)])+" +- "+str(consumpRatioErr[i][finalM[i].index(junk)])+'\n')
		# Sanity check
		if BR > MR:
		  print "ERROR: Maximal Respiration is lower than Basal Respiration!", conditions[i]
		# Calculate Spare Respiratory Capacity
		SRC = MR - BR
		SRCerr = numpy.sqrt(MRerr**2 + BRerr**2)
		gg.write('\t'+"Spare Respiratory Capacity = "+str(SRC)+" +- "+str(SRCerr)+'\n')
	else:
		print "ERROR: incorrect number of measurements. Only 14 and 12 are currently implemented."
		raise SystemExit

# Close output file
g.close()
