# RL-Xid Overview

This repository contains the code for all available versions of RL-Xid.
RL-Xid is an automated process which cross-identifies extended extragalactic 
radio sources with their optical counterparts using the innovative method of 
ridgelines.  The background and details of this process are described in my 
paper:

"The application of ridgelines in extended radio source cross-identification" 
(B. Barkus et al.) 2022, MNRAS 509, 1 found at: https://arxiv.org/abs/2110.05254

Thre are three main folders:

LOFAR
-----
This is the folder containing the main code running on the LOFAR Two

	DR1
The initial production and testing was run on the LOFAR DR1 datasets, and is 
documented in Barkus+22. This version of the code no longer needs to be run as it
is superceeded by the DR2 version.

	DR2
The current version is being implemented with LOFAR DR2, and is documented in
The LOFAR Two-metre Sky Survey (LoTSS) VI. Optical identifications for the 
second data release (M.J. Hardcastle et al., incl. B. Barkus) in prep. This is the
running version details in the readme.md for DR2.

MeerKAT
-------
There is a version adapted for use with a sample of MeerKAT MIGHTEE data that
is still under production. This does not have the DR2 improvements and is not currently
running.



