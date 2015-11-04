% Offest the stroke,
%   by subtracting offset from each time point
%
% Input
%  stk: [m x 2]
%  offset: [1 x 2]
%
% Output
%  stk: [m x 2]
function stk = offset_stk(stk,offset)
    assert(numel(offset)==2);
    n = size(stk,1);
    sub = repmat(offset,[n 1]);
    stk = stk - sub;
end