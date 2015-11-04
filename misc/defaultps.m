function PM = defaultps
    
    % Library to use    
    PM.libname = 'library';
    
    % Model parameters
    PM.K = 5; % number of particles to use in search algorithm

    % image model parameters
    PM.ink_ncon = 2; % number of convolutions
    PM.imsize = [105 105]; % image size
    PM.ink_pp = 2;       % amount of ink per point
    PM.ink_max_dist = 2; % distance between points to which you get full ink
    PM.ink_a = 0.5;      % ink parameter 1
    PM.ink_b = 6;        % ink parameter 2
    
    % Creating a trajectory from a spline
    PM.spline_max_neval = 200; % maxmium number of evaluations
    PM.spline_min_neval = 10;  % minimum
    PM.spline_grain = 1.5;     % 1 traj. point for every this many units pixel distance)
    
    % Max/min noise parameters for image model
    PM.max_blur_sigma = 16; % blur kernel width
    PM.min_blur_sigma = 0.5;
    PM.max_epsilon = 0.5;   % pixel flipping
    PM.min_epsilon = 1e-4;
    
    % search parameters
    PM.max_affine_scale_change = 2;  % scale changes must be less than a factor of 2
    PM.max_affine_shift_change = 50; % shift changes must less than this
    
    % MCMC PARAMETERS
    
    % details about the chain
    PM.mcmc.nsamp_type_chain = 200; % number of samples to take in the MCMC chain (for classif.)
    PM.mcmc.nsamp_type_store = 10; % number of samples to store from this chain (for classif.)
    PM.mcmc.nsamp_token_chain = 25; % for completion (we take last sample in this chain) 
    
    % mcmc proposal parameters (Note these are based on lib.tokenvar
    % parameters, although here they are hard-coded for convenience)
    PM.mcmc.prop_gpos_sd = 1; % global position move
    PM.mcmc.prop_shape_sd = 3 / 2; % shape move
    PM.mcmc.prop_scale_sd = 0.0235; % scale move
    PM.mcmc.prop_relmid_sd = 0.2168; % attach relation move
    PM.mcmc_prop_relpos_mlty = 2; % multiply the sd of the standard position noise by this to propose new positions from prior
    
end