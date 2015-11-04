# BPL model for one-shot learning

Bayesian Program Learning (BPL): Matlab source code for one-shot learning of handwritten characters.

### Citing this code
Please cite the following paper:


Lake, Brenden M., Salakhutdinov, Ruslan, and Tenenbaum, Joshua B. (in press). Human-level concept learning through probabilistic program induction. _Science_.


### Pre-requisites 

**The Lightspeed Matlab toolbox**   
(http://research.microsoft.com/en-us/um/people/minka/software/lightspeed/)   
Please install and add the lightspeed functions to your Matlab path.

**The Omniglot data set**   
(https://github.com/brendenlake/omniglot)   
Place these two Omniglot files in the 'data/' directory:   
matlab/data_background.mat   
matlab/data_evaluation.mat

### Using the code

**Setting your path**   
First, you must add all of the sub-directories to your Matlab path. While in the main BPL directory type this command:

```matlab
addpath(genpath('.'));
```

**Pre-processing stroke data**   
This only needs to be run once. From the 'data' directory, run:

```matlab
omniglot_preprocess;
```

This will create the 'data_background_processed.mat' and the 'data_evaluation_processed.mat' files for accessing the Omniglot dataset with pre-processed stroke data.

**Parsing demo**   
To run the parsing demo, type 

```matlab
test_parser;
```

### Computing resources

Most experiments will require a multi-CPU cluster to run in a reasonable amount of time. Fitting motor programs to images of characters can be run in parallel.


### Compatibility

Code was developed as tested on MATLAB R2013a and Lightspeed toolbox version 2.6.