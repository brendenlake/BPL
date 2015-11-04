%
% Get a motor trajectory from the b-spline control points,
% using an adaptive method to choose the number of evaluations
% based on the distance along the trajectory.
%
% Input
%
%  P: [ncpt x 2] control points
%  neval: (optional) number of evaluations
%    otherwise, we choose this adaptively
%
% Output
%  stk: [m x 2] trajectory
%
function stk = get_stk_from_bspline(P,neval)
    nland = size(P,1);
    if ~exist('neval','var')
        
        % set the number of evaluations adaptively,
        % based on the size of the stroke
        PM = defaultps;
        neval = PM.spline_min_neval;
        s = bspline_gen_s(nland,neval);
        stk = bspline_eval(s,P);
        sumdist = sum_pair_dist(stk);
        neval = max(neval,ceil(sumdist./PM.spline_grain));
        neval = min(neval,PM.spline_max_neval);
    
    end
    
    s = bspline_gen_s(nland,neval);
    stk = bspline_eval(s,P);
end

function s = sum_pair_dist(D)
    x1 = D(1:end-1,:);
    x2 = D(2:end,:);
    z = sqrt(sum((x1-x2).^2,2));
    s = sum(z);
end