Below are instructions to run RL-XiD on batches for a sky area within DR2. The main differences between the DR2 and DR1 versions of the code are to account for the different optical/IR catalogue inputs. The process below is designed so that the main steps are run in parallel via qsub multiple job submission on the UHHPC / LOFAR-UK cluster.

To run the full ridgeline ID process on a new sky area:

(1) make a subdirectory and put in it the relevant radio catalogue, components catalogue and optical catalogue. Components cat should be named "components-in.fits", and both radio and optical catalogue should have standard DR2 columns.
(2) Run rlcode/filter_batch.py with the radio and optical catalogue names as arguments. This filters the catalogues and then batches then appropriately for parallel processing.
(3) Run rlcode/opt_batch.py (at some point should incorporate this in filter_batch.py)
(3) Run DR2_setup.py in each batch subdirectory, using the run_setup.qsub file. Change FIL, DIR and python call within the qsub file to ensure correct files are referred to. This creates all the subdirectories needed and generates cutouts/npy arrays ready for input to RL-XiD.
(4) Run DR2_ridgelines.py on all the batches via an appropriate qsub file (e.g. run_ridgelines_fall.qsub) to generate the ridgelines.
(5) Run DR2_LR.py on all the batches via appropriate qsub fil (e.g. run_lr_fall.qsub) to carry out the likelihood ratio analysis
(6) Run unbatch.py to collate results and generate a test sample for quicklook checking. 

Step (6) outputs three files, a csv and fits hosts list, and a fits random subsample of 300 sources that can be used for visual checking. 

The csv allhosts file *has no LR threshold applied*, so includes the best host for all radio sources. It has columns:
Source_Name, RA, DEC, UNWISE_OBJID, UID_L from the input catalogues
LRMagBoth = likelihood ratio of best host
PossFail = flag for failed ridgelines

The fits 'allhosts' file is filtered to only include sources where the host has 'LRMagBoth>1' (this is the selected threshold for DR2 WEAVE samples, but could be modified). It doesn't include the columns UID_L and PossFail, but includes optRA_RLC and optDEC_RLC columns with the hosts RA and DEC (e.g. for use with visual inspection codes).

Finally the 'randhosts' fits file is a random subset of the 'allhosts' file containing 300 objects to enable easy visual checking. 
