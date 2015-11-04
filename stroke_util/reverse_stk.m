% Reverse the stroke
%
% Input
%  stk: [m x 2]
%
% Output
%  stk: [m x 2]
%
function stk = reverse_stk(stk)
    stk = stk(end:-1:1,:);
end