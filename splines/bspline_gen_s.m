% Generate time points for evaluating spline.
%
%  The convex-combination of the endpoitns with
%  five controls points are 80 percent the last cpt
%  and 20 percent the control point after that.
%
% Input
%   nland: number of landmarks
%   neval: number of evaluations
%
% Output
%   s: the time points used to evaluate spline
function [s,lb,ub] = bspline_gen_s(nland,neval)
    if nargin < 2
       neval = 200; 
    end
    
    lb = 2;
    ub = nland+1;
    epts = [lb ub];
    
    len = epts(2)-epts(1);
    int = len/(neval-1);
    s = epts(1):int:epts(2);    
end