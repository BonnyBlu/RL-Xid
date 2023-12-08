# MeerKAT Code Instructions

PLEASE DO NO TUSE THIS VERSION DOES NOT REALLY WORK! - GO LOOK AT DR2 IN LOFAR

This is a work in progress version of RL-Xid.

This version will draw ridgelines for data from the MIGHTEE survey, and has been tested on a small sample from the XMM-LSS dataset.

This version does not include the upgrades to LoTSS DR2 and therefore does not include:
- batching code for large datasets
- Healpix to divide the sky area in a sensible fashion (speeds the code up, and deals with edge cases which were causing the ID out  of region errors)
- Improved version of FloodFill/Masking function from J. Croston/B. Mingo (to remove the dependancy on PyRegion and speed the code up)
