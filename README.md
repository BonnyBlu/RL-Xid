# RL-Xid Overview

This repository contains the code for all available versions of RL-Xid. RL-Xid is an automated process which cross-identifies extended extragalactic radio sources with their optical counterparts using the innovative method of ridgelines.  The background and details of this process is described in my paper:

"The application of ridgelines in extended radio source cross-identification" 
(B. Barkus et al.) 2022, MNRAS 509, 1 found at: https://arxiv.org/abs/2110.05254

There are two main folders:

	LOFAR
This is the folder containing the main code running on the LOFAR Two-meter Sky Survey (LoTSS).

- DR1

The initial production and testing was run on the LoTSS DR1 datasets, and is documented in B. Barkus et al. 2022. This version of the code no longer needs to be run as it is superceeded by the DR2 version.

- DR2

The current version is being implemented with LoTSS DR2, and is documented in The LOFAR Two-metre Sky Survey (LoTSS) VI. Optical identifications for the second data release (M.J. Hardcastle et al., incl. B. Barkus) in prep. This is the current version and details are given in the README.md for DR2.

	MeerKAT

There is a version adapted for use with a sample of MeerKAT MIGHTEE data that is still under production. This does not have the DR2 improvements and is not currently running.


# Component Files

There are six component files that make up the main structure of RL-Xid.  These files are needed regardless of the telescope/survey the code is being run on. The process is split into two parts, and the ridgelines and likelihood ratio can be run sperately from each other, however the likelihood ratio requires the ridgelines to be drawn in order for it to run.

The files are listed and described in detail below, any additional files needed for specific surveys/telescopes are detailed in the appropriate README.md. Where a file is imported using an alias into another file the alias has been given as these are constistent across all versions of the code.  The placeholder XXX is used in the filenames to represent the location of version name indicator. For example, in XXX - LikelihoodRatio.ipynb the XXX can be replaced by DR1, DR2, or MK depending on if the file is being used for LoTSS DR1, LoTSS DR2, or MeerKAT.  This XXX is the part that will change when creating new versions and needs to be changed in the imports of each file.

The six files, and their import aliases are:
- XXX - LikelihoodRatio.ipynb
- XXX - Ridgelines.ipynb
- RLConstantsXXX.py 		imported as RLC
- RidgelineFilesXXX.py 		imported as RLF
- ridge_toolkitXXX.py
- SourceSearchXXX.py

The following is a description of each of these files including the contents of the file, which files are imported into each other, and the outputs they produce.

	XXX - LikelihoodRatio.ipynb

This is a jupyter notebook that runs the likelihood ratio part of RL-Xid. The ridgelines need to be drawn first before the likelihood ratios can be calculated.  The first cell will check to see if the ridgelines are present and run the code if needed. Each cell contains a description of what it is doing.

All imports are listed at the start, and this file requires RidgelineFilesXXX.py, ridge_toolkitXXX.py, and SourceSearchXXX.py to run.

The ouputs depend on which version is running and whether the likelihood ratios are being compared to a known dataset of hosts, or a threshold has been set to select the most probable counterpart.

	XXX - Ridgelines.ipynb

This is a jupyter notebook that runs the ridgeline drawing code of RL-Xid.  This is a short notebook containing the functions to select sources by flux and size from the given radio catalogue, create the cutouts to perform the ridgeline drawing process, and the ridgeline drawing function itself (this is a nested function of many steps).

All imports are listed at the start, and this file requires RidgelineFilesXXX.py, RLConstantsXXX.py, and ridge_toolkitXXX.py to run.

The outputs are the drawn ridgelines and the information text files, and a list of failed sources with images, see B. Barkus et al. 2022 for details.

	RLConstantsXXX.py

A list of the constants imported into the code files of RL-Xid. Created for ease of finding, changing, and adding constants to multiple large files.

	RidgelineFilesXXX.py

A list of the filenames used to open, save, and create files in RL-Xid. An overall path name has been created for naming ease, and to help build and maintain the correct file structure as many files will be produced in the process. This file also contains a list of the column names imported from the main catalogues, to allow for easy adjustment between muiltple different catalogues for alternate versions.

	ridge_toolkitXXX.py

This file contains all the functions required to run the ridgeline drawing and surface brightness processes.  Each function has a detailed doc string that describes what the function does, what parameters it takes, and what outputs it returns.

The imports are listed at the start, and this file requires RLConstantsXXX.py and RidglineFilesXXX.py to run.

	SourceSearchXXX.py

This file contains all the functions needed in the likelihood ratio process.  There is a description of the file at the start, and each function has a detailed doc string.  These describe what the function does, what parameters it takes, and what outputs it returns.

The imports are listed at the start, and this file requires RLConstantsXXX.py, RidgelineFilesXXX.py and ridge_toolkitXXX.py to run.


# Running RL-Xid

Currently the working version of RL-Xid to use is the one set up for LOFAR LoTSS DR2. The details on how to run this using qsub batching is given in the README.md. To run this without the batching, the batch code can be run one at a time, or the batch numbers can be set to produce a single group.
