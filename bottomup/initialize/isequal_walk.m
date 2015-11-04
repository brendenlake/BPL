%
% Are two random walks equal?
% This is a faster version of "isequal" 
%
function eq = isequal_walk(S1,S2)

    eq = true;
    
    % see if they have the same number of strokes
    if length(S1) ~= length(S2)
        eq = false;
        return
    end    
    
    for i=1:length(S1)
        
        % check same size, then check equality of the trajectory
        if ~isequal(size(S1{i}),size(S2{i})) || ~aeq(S1{i},S2{i})
            eq = false;
            return
        end
        
    end
    
end