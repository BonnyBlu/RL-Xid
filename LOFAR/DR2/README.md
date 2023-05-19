# DR2 Instructions

This is the current version of RL-Xid that is being implemented on LoTSS DR2. The files included here are the main six required to run RL-Xid and additional files to run everything in parallel batches. The following are the instructions for running RL-Xid in batch mode.  There are some modifications that need to be made to some of the batch files and these can be found in READMEFull.txt.  In order to run as batches DR2 - likelihoodRatios.ipynb and DR2 - Ridgelines.ipynb have been converted to .py scripts.

The original six files can still be run without the batching, however the batching files will provide:
- file structure creation
- faster code through better sky area selection
Therefore if a single batch is required it is better to adapt the batch numbers to provide a single batch rather than run the individual six files.

Below are instructions to run RL-XiD on batches for a sky area within
DR2. The main differences between the DR2 and DR1 versions of the code
are to account for the different optical/IR catalogue inputs. The
process below is designed so that the main steps are run in parallel
via qsub multiple job submission on the UHHPC / LOFAR-UK cluster.

To run the full ridgeline ID process on a new sky area:

1. make a subdirectory and put in it the relevant radio catalogue,
   components catalogue and optical catalogue. Components cat should
   be named `components.fits`, and both radio and optical catalogue
   should have standard DR2 columns.
2. Run `healpix_batch.py` with the radio catalogue, optical catalogue
   and component catalogue as arguments. This filters the catalogues
   and then batches then appropriately for parallel processing. It
   will generate `N` directories with names `hp_0`, `hp_1`... 
3. Run `DR2_setup.py` in each batch subdirectory, using the
   `run_setup_hp.qsub` file. Use the `-v RDIR` option to specify the
   root directory and the `-t` option to specify the range up to the
   largest numbered `hp` directory. This creates all the subdirectories
   needed and generates cutouts/npy arrays ready for input to RL-XiD.
4. Run `DR2_ridgelines.py` on all the batches via an appropriate qsub
   file (e.g. `run_ridgelines_hp.qsub`) to generate the ridgelines.
5. Run `DR2_LR.py` on all the batches via appropriate qsub file
   (e.g. `run_lr_hp.qsub`) to carry out the likelihood ratio analysis
6. Run `unbatch.py` to collate results and generate a test sample for
   quicklook checking.

Step 6 outputs three files, a csv and fits hosts list, and a fits
random subsample of 300 sources that can be used for visual checking.

The `allhosts.csv` file *has no LR threshold applied*, so includes the
best host for all radio sources. It has columns: `Source_Name`, `RA`,
`DEC`, `UNWISE_OBJID`, `UID_L` from the input catalogues, `LRMagBoth` =
likelihood ratio of best host, `PossFail` = flag for failed ridgelines.

The `allhosts.fits` file is filtered to only include sources where the host has 'LRMagBoth>1' (this is the selected threshold for DR2 WEAVE samples, but could be modified). It doesn't include the columns `UID_L` and `PossFail`, but includes `optRA_RLC` and `optDEC_RLC` columns with the hosts' `RA` and `DEC` (e.g. for use with visual inspection codes).

Finally the `randhosts.fits` file is a random subset of the `allhosts` file containing 300 objects to enable easy visual checking. 

This depends on the `dr2_catalogue` repository and `sizeflux_tools.py` must be on your Pythonpath.

Example run (assumes relevant packages in your `git` directory):

```
python ~/git/lotss-catalogue/dr2_catalogue/make_ridgeline_input.py
mkdir ridgeline
mv ridgeline-in.fits ridgeline
cd ridgeline
ln -s ../components-v0.6.fits components.fits
python ~/git/RL-Xid/LOFAR/DR2/healpix_batch.py ridgeline-in.fits ../optical.fits components.fits

qsub -t 0-366 -v RDIR=/beegfs/lofar/mjh/rgz/Spring/ridgeline ~/git/RL-Xid/LOFAR/DR2/example\ qsub\ files/run_setup_hp.qsub

qsub -t 0-366 -v RDIR=/beegfs/lofar/mjh/rgz/Spring/ridgeline ~/git/RL-Xid/LOFAR/DR2/example\ qsub\ files/run_ridgelines_hp.qsub

qsub -t 0-366 -v RDIR=/beegfs/lofar/mjh/rgz/Spring/ridgeline ~/git/RL-Xid/LOFAR/DR2/example\ qsub\ files/run_lr_hp.qsub

python3.6 ~/git/RL-Xid/LOFAR/DR2/unbatch.py hp_ 
```
