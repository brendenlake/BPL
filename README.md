# BPL model for one-shot learning

Matlab source code for one-shot learning of handwritten characters with Bayesian Program Learning (BPL).

### Citing this code
Please cite the following paper:

[Lake, B. M., Salakhutdinov, R., and Tenenbaum, J. B. (2015). Human-level concept learning through probabilistic program induction.](http://www.sciencemag.org/content/350/6266/1332.short) _Science_, 350(6266), 1332-1338.


### Pre-requisites 

**Matlab Toolboxes**   
Optimization Toolbox   
Statistics Toolbox (before R2015a) OR Statistics and Machine Learning Toolbox   
Image Processing Toolbox   
Curve Fitting Toolbox

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
addpath(genpath(pwd));
```

**Pre-processing stroke data**   
This only needs to be run once, and it can take up to 5 minutes to complete. From the 'data' directory, run:

```matlab
omniglot_preprocess;
```

This will create the 'data_background_processed.mat' and the 'data_evaluation_processed.mat' files for accessing the Omniglot dataset with pre-processed stroke data.

**Parsing demo**   
To run the model fitting demo, type 

```matlab
demo_fit;
```

**One-shot classification**   
First, download the pre-computed models and unzip such that 'model_fits' and 'model_refits' are sub-directories of the 'classification' directory.   
http://cims.nyu.edu/~brenden/supplemental/BPL_precomputed/model_fits.zip   
http://cims.nyu.edu/~brenden/supplemental/BPL_precomputed/model_refits.zip

To run the model re-fitting demo, enter the 'classification' directory and type:

```matlab
demo_refit;
```

To measure classification error rate with the pre-computed results, type:

```matlab
run_classification;
```

**One-shot exemplar generation**

First, enter the 'generate_exemplars' directory and unzip 'model_fits.zip' so 'model_fits' is now a sub-directory. These are the pre-computed models that were used in the exemplar generation visual Turing tests.

To run the exemplar generation demo, from within the 'generate_exemplars' directory type:

```matlab
demo_generate_exemplar;
```

### Computing resources

Most experiments will require a multi-CPU cluster to run in a reasonable amount of time. Fitting motor programs to images of characters can be run in parallel.

The parsing and one-shot classification demos include a 'fast_mode' option (on by default) which allows for the demo to run quickly, skipping the expensive procedure of fitting the strokes to the details of the image. Use this mode with caution, as this mode is much cruder and was not used in the paper results.

### Compatibility

Code was developed and tested on MATLAB R2013a and Lightspeed toolbox version 2.6.