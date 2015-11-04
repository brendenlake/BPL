% Get the center of mass of a stroke
%
% Input
%  stk: [m x 2]
%
% Output
%  com: [scalar]
%
function com = com_stk(stk)
    com = mean(stk,1);
end