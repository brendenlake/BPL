% Fit a motor program to an image.

K = 1; % number of parses we want to collect
verbose = true; % describe progress and visualize parse
load('L_img','img');
G = parser(img,K,verbose);