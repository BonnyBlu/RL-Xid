# RL-XID Overview

This repository contains the code for all available versions of RL-Xid. RL-Xid is an automated process which cross-identifies extended extragalactic radio sources with their optical counterparts using the innovative method of ridgelines.  The background and details of this process are described in my paper:

"The application of ridgelines in extended radio source cross-identification" 
(B. Barkus et al.) 2022, MNRAS 509, 1 found at: https://arxiv.org/abs/2110.05254

There are two main folders:

	LOFAR
This is the folder containing the main code running on the LOFAR Two-meter Sky Survey (LoTSS).

- DR1

The initial production and testing was run on the LOFAR DR1 datasets, and is documented in Barkus+22. This version of the code no longer needs to be run as it is superceeded by the DR2 version.

- DR2

The current version is being implemented with LOFAR DR2, and is documented in The LOFAR Two-metre Sky Survey (LoTSS) VI. Optical identifications for the second data release (M.J. Hardcastle et al., incl. B. Barkus) in prep. This is the current version and details are given in the readme.md for DR2.

	MeerKAT

There is a version adapted for use with a sample of MeerKAT MIGHTEE data that is still under production. This does not have the DR2 improvements and is not currently running.


# Component Files

There are six component files that make up the main structure of RL-Xid.  These files are needed regardless of the telescope/survey the code is being run on. The process is split into two parts and the ridgelines and likelihood ratio can be run sperately from each other, however the likelihood ratio requires the ridgelines to be drawn in order for it to run.

The files are listed and described in detail below, any additional files needed for specific surveys/telescopes are detailed in the appropriate readme.md. Where a file is imported using an alias into another file the alias has been given as these are constistent across all versions of the code.  The placeholder XXX is used in the filenames to represent the location of version name indicator. For example, in XXX - LikelihoodRatio.ipynb the XXX can be replaced by DR1, DR2, or MK depending on if the file is being used for LoTSS DR1, LoTSS DR2, or MeerKAT.

The six files, and their import aliases are:
- XXX - LikelihoodRatio.ipynb
- XXX - Ridgelines.ipynb
- RLConstantsXXX.py 		imported as RLC
- RidgelineFilesXXX.py 		imported as RLF
- ridge_toolkitXXX.py
- SourceSearchXXX.py




# Running RL-Xid
