% Get the overall center-of-mass for
%  a character
% 
% Input
%   char: [ns x 1 cell] array of strokes
% 
% Output
%   COM: [scalar] center of mass
function COM = com_char(char)

    ns = length(char);
    lens = zeros(ns,1);
    wsum = zeros(ns,2);
    for i=1:ns
       stk = char{i};
       lens(i) = length_stk(stk);
       wsum(i,:) = com_stk(stk) .* lens(i);
    end
    
    COM = sum(wsum,1) ./ sum(lens);
    assert(numel(COM)==2);
end