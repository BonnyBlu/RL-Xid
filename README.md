# RL-Xid Overview

This repository contains the code for all available versions of RL-Xid.
RL-Xid is an automated process which cross-identifies extended extragalactic 
radio sources with their optical counterparts using the innovative method of 
ridgelines.  The background and details of this process are described in my 
paper:

"The application of ridgelines in extended radio source cross-identification" 
(B. Barkus et al.) 2022, MNRAS 509, 1 found at: https://arxiv.org/abs/2110.05254

- The initial production and testing was run on the LOFAR DR1 datasets, and is 
documented in Barkus+22.
- The current version is being implemented with LOFAR DR2, and is documented in
The LOFAR Two-metre Sky Survey (LoTSS) VI. Optical identifications for the 
second data release (M.J. Hardcastle et al., incl. B. Barkus) in prep.
- There is a version adapted for use with a sample of MeerKAT MIGHTEE data that
is still under production.

The structure of this repository is as follows, and each folder contains a detailed readme:
RL-Xid ---> LOFAR ---> DR1: Testing code
	|     |
	|     -------> DR2: Current code
	|
	---> MeerKAT: In progress code
