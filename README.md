SCR Model Code Version 0.1

This document:
  I. Files
  II. How to run the code
  III. Limitations
  IV. How to test the code
  V. Release information

I. Files
------------
1. modelSCR_core.py contains the core modeling routines that implement the model described in the SCR SC 2010 paper (see References section below).

2. modelSCR.py is a driver program that sets up multi-parameter experiments based on command line arguments. For each experiment specified, it calls the routines in modelSCR_core.py.


II. How to run the code
--------------------------

For a full parameter list, run:
$  modelSCR.py -h

The current driver program, modelSCR.py, takes parameters that tell it what experiments to run, based on inputs of: checkpoint interval time, machine to emulate, overhead factor, failure factor, and L1 and L2 checkpoint counts.
 - checkpoint interval time: This is the time between the ending of the last checkpoint write and the beginning of the next checkpoint write. The driver program takes a low value (--checkIntervalLow), a high value (--checkIntervalHigh) and a step value (--checkIntervalStep) and will run an experiment for each of these values for each of the other specified parameters.
 - machine to emulate: There are values hard-coded into the driver program for the failure rates and overheads of real machines. By specifying a machine name, you specify what failure and overhead data to use (--machine). Current supported machines are coastal, atlas, and "made_up" (not real data).
 - overhead factor: Specify a value to multiply the highest level checkpoint overhead by, e.g. if the specified overhead factor is 2 and the highest checkpoint level is Lustre, then the time to take and recover from Lustre checkpoints will be multiplied by 2. (-o, --overheadFactor)
 - failure factor: Specify a value to multiply the failure rates by, e.g. if the specified failure factor is 2, then all the failure rates will be multiplied by 2 (-f, --failureFactor)
 - L1 and L2 checkpoint counts: Specify how many level 1 and level 2 checkpoints are taken. These parameters take low values (L1Low, L2Low), high values (L1high, L2high), and step values (L1Step,L2Step). The driver runs experiments for each combination of L1 and L2 values for each of the other parameters specified.


--expDir         Specify the name of the directory to put the data in. Summary files (see directly below for discussion of summary files) will be named: maxes.kX.oY.fZ.data, where X is the number of levels in the system, Y is the overhead factor, and Z is the failure factor. Full results (see directly below for discussion of full results) will be in files called: A_intervalTime_B_C_D_E_F_G_H_I_J.data, where 
A is the expPrefix (see below) 
B is the low checkpoint interval value
C is the high checkpoint interval value
D is the checkpoint interval step
E is the number of checkpoint levels
F is a dummy value
G is the level 1 count
H is the level2 count
I is the overhead factor
J is the failure factor.

Optional arguments:

--fullResults      By default, it will only write summary files of the maximum efficiencies for each of the parameter combinations over the span of checkpoint interval times, i.e. for each l1,l2,o,f it finds the checkpoint interval that results in the highest efficiency. If you specify --fullResults, it will print the results of each experiment for each time interval.
--epsilon         To account for floating point round-off error, when comparing values with 0, we use an epsilon instead. You can specify a different value
--debug           Turn on debug checks and printing
--expPrefix       Specify a prefix to append to data files

Example: 

$   ./modelSCR.py --L1Low 1 --L1High 3 --L1Step 1 --checkIntervalLow 1 --checkIntervalHigh 10 --checkIntervalStep 1 -f 1 -o 1 --machine atlas --expDir testing

This will simulate running SCR on the machine atlas and will write the data files to the directory "testing". It will run with level 1 checkpoint counts of 1,2,3 and checkpoint intervals of 1, 2, ..., 10.  It will only write summary files of the maximum efficiencies.


III. Limitations
------------------
1. We only have real data for machines that take checkpoints at 2 levels. Therefore the L2 (level 2) command line arguments are ignored at this time. The "made-up" machine has 3 levels of failure and overhead inputs. These values are hard-coded in, and the level counts are altered accordingly. A better solution would be to use input decks or similar.

2. The code is hard-coded to only work for systems with 2 or 3 levels. To do generic experiments (1 <= l <= L), the driver code will have to be altered.


IV. How to test the code
--------------------------
If you make changes, you can test the code by running regression tests. The script regressTest.py will run a series of tests, compare the results with what are believed to be good results, and then print out success or failure messages.


V. Release
-------------
LLNL-CODE-817266
See the LICENSE file for license information

VI. References
----------------
Adam Moody, Greg Bronevetsky, Kathryn Mohror, Bronis R. de Supinski, "Design, Modeling, and Evaluation of a Scalable Multi-level Checkpointing System," LLNL-CONF-427742, Supercomputing 2010, New Orleans, LA, November 2010. https://dl.acm.org/doi/10.1109/SC.2010.18

Kathryn Mohror, Adam Moody, Greg Bronevetsky, Bronis R. de Supinski, "Detailed Modeling and Evaluation of a Scalable Multilevel Checkpointing System," in Transactions on Parallel and Distributed Systems, LLNL-JRNL-564721, Transactions on Parallel and Distributed Systems, 25(9):2255-2263, Sept. 2014. https://ieeexplore.ieee.org/document/6494566
