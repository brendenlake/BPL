% Rescale a stroke,
%   by multiplying the x and y-dimension by 
%   scale(1) and scale(2), respectively
%
% Input
%  stk: [m x 2]
%  scale: [1 x 2]
%
% Output
%  stk: [m x 2]
%
function stk = rescale_stk(stk,scale)
    assert(numel(scale)==2);
    stk(:,1) = stk(:,1) .* scale(1);
    stk(:,2) = stk(:,2) .* scale(2);
end