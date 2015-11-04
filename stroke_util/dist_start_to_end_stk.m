% Get distance of a straight line from the start to end 
% point of a stroke.
% 
% Input
%  stk: [m x 2]
% 
% Output
%   dist: [scalar] distance from starting position to ending position
function dist = dist_start_to_end_stk(stk)
    mystart = stk(1,:);
    myend = stk(end,:);
    dist = norm(mystart-myend);
end