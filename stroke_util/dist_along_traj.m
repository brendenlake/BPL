% Get sum distance along the trajectory, from beginning to end
% 
% Input
%  stk: [m x 2]
% 
% Output
%   dist: [scalar] distance from starting position to ending position
function dist = dist_along_traj(stk)
    [m,dim] = size(stk);
    assert(dim==2);
    dist = 0;
    for i=2:m
       dist = dist + norm(stk(i,:)-stk(i-1,:));        
    end
end