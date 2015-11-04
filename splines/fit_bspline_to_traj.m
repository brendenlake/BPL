%
% Fit a b-spline to "stk" with "nland" landmarks
%
function P = fit_bspline_to_traj(stk,nland)
    neval = length_stk(stk);
    s = bspline_gen_s(nland,neval);
    P = bspline_fit(s,stk,nland);   
end