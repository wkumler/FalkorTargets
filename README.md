# FalkorTargets

A repo to share with Katherine to review my annotation of targeted compounds in the lab

## Contents

### Scripts

**stan_assignments.R**: The main script, used to annotate features with compound names according to rules
 - Requires real_peaks.csv, addiso_peaks.csv, and falkor_stans.csv
 - Outputs stan_assignments.csv, connects feature numbers to compound names

**stan_vis.Rmd**: An Rmarkdown document used to visualize the peaks and manually check their annotations
 - Requires stan_assignments.csv and stan_data.rds
 - Outputs stan_vis.html, a knit document with plotted EICs and metadata.

**stan_data_maker.R**: Creates stan_data.rds by opening each MS file, extracting relevant data, and formatting. 
Shouldn't need to be run and will break if run not on Will's computer. Fix the file paths to point to the 
correct files if you'd like to use this

### Data

**real_peaks.csv**: The list of peaks found by peakpicking, excluding adducts and isotopes

**addiso_peaks.csv**: The list of peaks found by peakpicking to be adducts or isotopes

**falkor_stans.csv**: A clean version of the Ingalls Lab Standards sheet

**stan_assignments.csv**: Data frame associating each standard with a feature number as found by peakpicking

**stan_data.rds**: A data table with EIC information (mz, rt, int) for each standard to be plotted

### Output

**stan_vis.html**: A knit document with plotted peaks and associated metadata

## Regina comments ##
I like your layout of Scripts, data, output. I might clarify whether things need to be run in a particular order, or if they can just be run independently than I would mention that too.
