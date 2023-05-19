# MeerKAT Code Instructions

This is a work in progress version of RL-Xid.

This version will draw ridgelines for data from the MIGHTEE survey, and has been tested on a small sample from the XMM-LSS dataset.

This version does not include the upgrades to LoTSS DR2 and therefore does not include:
- batching code for large datasets
- Healpix to divide the sky area in a senisble fashion (speed the code up and deal with edge cases causing ID out region errors)
- Improved version of FloodFIll/Mask function from J. Croston/B. Mingo (to remove the dependancy on PyRegion and speed the code up)
