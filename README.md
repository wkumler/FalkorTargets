# FalkorTargets

A repo to share with Katherine to review my annotation of targeted compounds in the lab

## Contents

### Scripts

stan_assignments.R: The main script, used to annotate features with compound names according to rules
 - Requires real_peaks.csv, addiso_peaks.csv, and falkor_stans.csv
 - Outputs stan_assignments.csv, a cleaned version of falkor_stans.csv annotated with feature numbers

stan_vis.Rmd: An Rmarkdown document used to visualize the peaks and manually check their annotations
 - Requires stan_assignments.csv and stans_data.rds
 - Outputs stan_vis.html, a knit document with plotted EICs and metadata.

stans_data_maker.R: Creates stans_data.rds by opening each MS file, extracting relevant data, and formatting. 
Shouldn't need to be run and will break if run not on Will's computer. Fix the file paths to point to the 
correct files if you'd like to use this

### Data

real_peaks.csv: The list of peaks found by peakpicking, excluding adducts and isotopes

addiso_peaks.csv: The list of peaks found by peakpicking to be adducts or isotopes

falkor_stans.csv: A clean version of the Ingalls Lab Standards sheet

stan_assignments.csv: Essentially a clean version of falkor_stans.csv, where each standard is 
annotated with a feature number as found by peakpicking

stans_data.rds: A data table with EIC information (mz, rt, int) for each standard to be plotted