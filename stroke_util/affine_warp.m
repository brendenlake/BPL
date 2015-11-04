%
% Affine warp defined by
% A(1) * [x y] + [A(2) A(3)]
% OR
% A(1:2) .* [x y] + [A(3) A(4)]
%
% Input
%  stk [n x 2] stroke
%  A : [3x1 or 4x1] affine warp
%
function stk = affine_warp(stk,A)
    
    n = size(stk,1);
    if numel(A)==3
        stk = A(1).*stk;
        stk = stk + repmat(A(2:3)',[n 1]);
    elseif (numel(A)==4)
        stk = stk .* repmat(A(1:2)',[n 1]);
        stk = stk  + repmat(A(3:4)',[n 1]);
    else
        error('invalid affine warp');
    end
end