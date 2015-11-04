% Get the minimum values in a stroke
%
% Input
%  stk: [m x 2]
%
% Output
%  mn: [1 x 2]
%
function mn = min_stk(stk)
    mn = min(stk,[],1);
    assert(numel(mn)==2);
end