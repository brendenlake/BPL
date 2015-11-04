%
% Parameters for the bottom-up parsing algorithm.
%
function ps = defaultps_bottomup
    
    ps = struct;
    
    % FOR THE RANDOM WALK MODEL
    ps.max_nwalk = 150; % maximum number of random walks we are willing to generate
    ps.max_nstroke = 100; % maximum number of strokes we are willing to parse
    
    ps.nwalk_det = 5; % how many times should we try a deterministic walk?
    ps.nwalk_order_mult = 2; % compared to "ninit", the number of parses we want,
                             % how many of the best ones should we optimize
                             % the stroke order of?
    ps.max_len = 8; % maximum number of strokes to allow
    ps.rad_junction = 5; % pixel distance radius around junction for computing angle

    % Discrete distribution on randomness parameters, which control
    % the stochastic nature of the search
    ps.int_exp_wt = [2 4 6]; % value determines how strongly we start at nodes with low degree
    ps.int_lambda_soft = [2 4 8 16]; % value determines how strongly we follow the closest angle
    
    % faux angles: moves of this type have this value for their "angle"
    ps.faux_angle_lift = 45; % value of lifting one's pen
    ps.faux_angle_repeat = 90; % value of repeating a previous edge
    ps.faux_angle_too_short = 45; %
    
    % FOR SPLINE SMOOTHING
    ps.int_rho = [1e-4,1e-2,1e-1,1]; % values to try for smoothing strokes,
                                     % starting with a straight line

    % CONVERTING TRAJECTORIES TO SUB-STROKES
    ps.traj_abs_error_lim = 3; % maximum acceptable error along spline fit for any point
    ps.ntry_split = 7; % try this many splits at each iteration of sub-stroke parsing
    ps.sigma_wiggle = 2; % s.d. of number of points to wiggle splits
     
    % FINAL CONVERSION TO MOTORPROGRAM
    ps.max_ns_all_perm = 5; % when searching for optimal stroke order, 
                            % try all permutations up to this number of strokes 
    ps.init_epsilon = .01;  % initial value of epsilon noise parameter
    ps.init_blur_sigma = 1; % initial value of the blur parameter

end