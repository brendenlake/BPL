% Get the overall range (x,y) that a character spans
% 
% Input
%   char: [ns x 1 cell] array of strokes
% 
% Output
%   range: [1 x 2 scalar] range of each dimension
function range = range_char(char)

    ns = length(char);
    mymax = zeros(ns,2);
    mymin = zeros(ns,2);
    
    % for each stroke
    for i=1:ns
       stk = char{i};
       mymax(i,:) = max_stk(stk);
       mymin(i,:) = min_stk(stk);
    end
    
    % Get overall max and min
    mx = max(mymax,[],1);
    mn = min(mymin,[],1);

    % Compute the range
    range = mx - mn;
    assert(numel(range)==2);
end