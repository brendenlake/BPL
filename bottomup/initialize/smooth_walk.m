%
% Smooth the strokes in a random walk parse
% of the image.
%
% Approach: Start with the smoothest possible stroke,
% and if that doesn't work, keep adding complexity
% 
% Input
%   S: [ns x 1 cell] random walk in IMAGE space
%   I : binary image
%   bool_viz : 
%
function S = smooth_walk(S,I,bool_viz)

    if ~exist('bool_viz','var')
       bool_viz = false; 
    end
    assert(iscell(S));
    
    ps = defaultps_bottomup;
    ns = numel(S);
    for sid=1:ns
       S{sid} = fit_smooth_stk(S{sid},I,bool_viz,ps);
    end
    
end