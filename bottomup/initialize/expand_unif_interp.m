%
% Expand a trajectory to achieve uniform spacing of "desired_dist"
% between all points. Only approximately achieves this spacing.
%
% Input
%   traj: [n x 2] stroke trajectory
%   desired_dist: [scalar] how far apart we want the items to be
%
% Output
%   build_traj: [m x 2] modified trajectory
%
function build_traj = expand_unif_interp(traj,desired_dist)
    
    [n,dim] = size(traj);
    assert(dim==2);
    
    % compute all pair-wise distances
    d = zeros(n-1,1);
    for i=2:n
        d(i) = norm(traj(i,:)-traj(i-1,:));
    end
    
    newsteps = cell(n,1); % extra steps
    for i=2:n
        
        int_divide = get_opt_spacing(d(i),desired_dist);
        
        % linear interpolation
        ntot = int_divide+1;        
        xi = linspace(0,1,ntot);        
        Y = [traj(i-1,:); traj(i,:)];        
        yi = interp1([0 1],Y,xi);
        
        % extract only the new steps
        newsteps{i} = yi(2:end-1,:);
    end
    
    % Build the trajectory with the new intervening segments
    build_traj = [];
    for i=1:n
        build_traj = [build_traj; newsteps{i}; traj(i,:)];
    end


end

% current distance d vs. optimal distance desired_dsit
% 
% Output
%  int_divide: how many sub-segments we want in this segment
function int_divide = get_opt_spacing(d,desired_dist)

    assert(isscalar(d));
    
    % compute candidate expansion factors, or
    % how many sub-segments we want in this length of traj.
    rough_int = d./desired_dist;
    floor_int = floor(rough_int);
    ceil_int = ceil(rough_int);
    
    % make sure we don't shrink the trajectories
    floor_int = max(floor_int,1);
    ceil_int = max(ceil_int,1);
    
    % compute error of those expansion factors,
    % in terms of missed distance for each inteval
    floor_error = floor_int*( abs(desired_dist-d/floor_int) );
    ceil_error = ceil_int*( abs(desired_dist-d/ceil_int) );
    
    if floor_error <= ceil_error
        int_divide = floor_int;
    else
        int_divide = ceil_int; 
    end
    
end