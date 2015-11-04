% Fit a uniform, cubic B-spline
%
% Input
%   sval: vector [k x 1] where 0 <= sval(i) <= n
%   cpts: [n x 2] array of control points
%
% Output
%    y: vector [k x 2] which is output of spline
function [y,Cof] = bspline_eval(sval,cpts)
    assert(isvector(sval));
    sval = sval(:);
    L = size(cpts,1);
    ns = length(sval);
    y = zeros(ns,2);
    
    S = repmat(sval,[1 L]);
    I = repmat(0:L-1,[ns 1]);
    Cof =  vectorized_bspline_coeff(I,S);
    sumC = sum(Cof,2);
    Cof = Cof ./ repmat(sumC,[1 L]);
    y(:,1) = Cof*cpts(:,1);
    y(:,2) = Cof*cpts(:,2);
    
    % Error checking
    assert(~any(isnan(y(:))));
end