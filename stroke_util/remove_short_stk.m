%
% Remove a subset of the strokes that are too short
%
% Input
%  nested: nested cell array of strokes, where the 
%       bottom level are arrays for each stroke
%  minlen: minimum length required for removal
%  mindist: minimum distance required for removal
% 
% Output
%  nested: nested cell array with elements removed
%
function nested = remove_short_stk(nested,minlen,mindist)   
    nested = process(nested,minlen,mindist);
end

function N = process(nested,minlen,mindist)

    m = numel(nested);
    
    % If we are at the character level
    if ~iscell(nested{1})    
        
        len = zeros(m,1);
        dist = zeros(m,1);
        for i=1:m
           stk = nested{i};
           len(i) = length_stk(stk);
           dist(i) = dist_along_traj(stk);
        end
        tormv = (len < minlen & dist < mindist);
        
        N = nested;
        if true %~all(tormv)
            N(tormv) = []; % clear small strokes
        end
        
        return
    end
    
    % Recursive call
    N = nested;
    for i=1:m
       N{i} = process(nested{i},minlen,mindist); 
    end
    
end