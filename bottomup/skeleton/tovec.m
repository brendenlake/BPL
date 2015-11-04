% Convert a cell array of vectors
% into a single long concat.'d vector
%
% Ouput
%  v [n x 1]: vector of all the cell array vectors combined
%  indx [n x 1]: which cell array each element came from
function [v,indx] = tovec(cell_of_vec)
    v = [];
    indx = [];
    for i=1:numel(cell_of_vec)
       els = cell_of_vec{i};
       if ~isempty(els)
           v = [v; vec(els)];
           indx = [indx; i*ones(length(els),1)];
       end
    end
end