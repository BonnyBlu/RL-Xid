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




Instructions
------------

The following are detailed instructions on how to adapt and run the files.  There are pathways that need to be changed and column names that need to be checked for your catalogues.  This is a comprehensive list of those checks that need to take place.


File structure
--------------

This structure will allow for multiple runs of the code with the results from each going into a new folder.  As long as the correct folder names are entered as recommended below.

RL-Xid		--------------- 'Code (folder)' - name as wish
(Main Folder)		|
			-------	'Results1 (folder)' - name as wish
			|
			-------	'Results2 (folder)' - etc


filter_batch.py
---------------

- line 15 - need to change bpwd to the pathway to the results folder.
- Check column names:
   - line 22 - the total flux of the source
   - line 23 - the width or major axis of the source
   - line 24 - the size or LGZ size of the source
- line 40 - change the number of sources per batch
- takes pathway to radio catalogue, pathway to optical catalogue, 'batch(folder)' as the input (create a batch name for the results)


**	run:	python filter_batch.py 'radio cat pwd' 'optical cat pwd' 'batch(folder)'	**


opt_batch.py
------------

- Check column names:
    - lines 24, 25, 32 and 33 - RA from the opt catalogue
    - lines 26, 27, 32 and 33 - Dec from the opt catalogue
- takes ../'results(folder)'/'batch(folder)' as input


** 	run:	python opt_batch.py ../'results(folder)'/'batch(folder)' 	**


RLConstantsDR2.py
-----------------

- Only update if you need too


RidgelineFilesDR2.py
--------------------

- line 17 - update pathway to component catalogue
- lines 32-37 - update any filenames to personal choice
- Check column names:
    - lines 55-61 - update to match radio catalogue (typically SASS)
    - lines 62-69 - update to match component catalogue (typically CSN and CCN)
    - lines 70-80 - update to match optical catalogue (typically OptMagA, OptMagP, PossRA, PossDEC)
- line 94 - change pathway to HostMagnitude_Info.txt
- line 95 - change pathway to HostMagnitude_InfoFull.txt


DR2_setup.py
------------

- Check column names:
    - line 79 - the width or major axis of the source
    - line 80 - the size or LGZ size of the source
    - lines 81, 83, 84, 85, 88 and 89 - update to match radio catalogue
    - lines 133, 134, 137 and 138 - update to match component catalogue


run_setup.qsub
--------------

- line 18 - change pathway to location of image mosaic
- line 19 - change pathway to version of python containing LOFAR packages
- line 20 - change pathway to the 'batch(folder)'${PBS_ARRAYID}
- line 23 - change pathway to the 'batch(folder)'$[PBS_ARRAYID}
- line 28 - change first pathway to include correct 'code(folder)'.  Change second pathway to the component catalogue

**	run: 	qsub -t 0-#(number of batches) run_setup.qsub	**


run_ridgelines.qsub
-------------------

- line 18 - change to correctly load python 3
- line 20 - change pathway to 'code(folder)'
- line 22 - change pathway to 'batch(folder)'${PBS_ARRAYID}

**	run:	qsub -t 0-#(number of batches) run_ridgelines.qsub		**


run_lr_fall.qsub
----------------

- line 19 - change to correctly load python 3
- line 21 - change pathway to 'code(folder)'
- line 23 - change pathway to 'batch(folder'${PBS_ARRAYID}

**	run:	qsub -t 0-#(number of batches) run_lr_fall.qsub	**


unbatch.py
----------

- line 28 - change number of random sample to chosen amount
- line 31 - change threshold for selection if desired
- line 32 - change number of random sample to chosen amount
- takes ../'results(folder)'/'batch(folder)' as input

**	run:	python unbatch.py ../'results(folder)'/'batch(folder)'	**




