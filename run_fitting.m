
% Fit a motor program to an image.
K = 3; % number of parses we want to collect
verbose = true; % describe progress and visualize parse
include_mcmc = false;
fast_mode = true;

load('I_img','img');
% load('Phi_img','img');
G = fit_motorprograms(img,K,verbose,include_mcmc,fast_mode);