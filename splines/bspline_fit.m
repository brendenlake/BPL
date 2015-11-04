% Fit a bspline using least-squares
% 
% Input
%   sval: time points (N x 1)
%   X: data points (N x 2)
%   L: number of control points to fit
% 
% Output:
%   P: (L x 2) optimal control points
function P = bspline_fit(sval,X,L)
    sval = sval(:);
    ns = length(sval);
    assert(isequal(size(X),[ns 2]));
   
    S = repmat(sval,[1 L]);
    I = repmat(0:L-1,[ns 1]);
    A = vectorized_bspline_coeff(I,S);
    
    sumA = sum(A,2);
    Cof = A ./ repmat(sumA,[1 L]);
    P = (Cof'*Cof)\Cof'*X;
    %P = inv(Cof'*Cof)*Cof'*X;
end