%
% Parameters for data pre-processing steps
%
function ps = defaultps_preprocess

    ps = struct;
    ps.tint = 50; % interpolate stroke so that we have datapoints every
                  % interval of this many milliseconds
                  
    ps.dthresh = 1; % if this much distance (norm) is not covered
                    % at each time point, then it's a pause
                    
    ps.max_sequence = inf; % maximum length of a stop sequence, before it is
                           % called a different sub-stroke
                      
    ps.space_int = 1; % we want approximately this much distance (norm)
                      % covered between successive points for spatial
                      % interpolation
end