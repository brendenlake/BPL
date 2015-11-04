% Get the maximum values in a stroke
%
% Input
%  stk: [m x 2]
%
% Output
%  mx: [1 x 2]
%
function mx = max_stk(stk)
    mx = max(stk,[],1);
    assert(numel(mx)==2);
end