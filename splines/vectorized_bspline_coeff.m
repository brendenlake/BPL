% function vectorized_bspline_coeff
% ------
% 
% Input
%  vi [n x m] 
%  vs: [n x m]
% 
%  see Kristin Branson's "A Practical Review of Uniform B-splines"
%
% Output
%   C [n x 1]: the coefficients
function C = vectorized_bspline_coeff(vi,vs)
    
    assert(isequal(size(vi),size(vs)));
    
    % Go through conditions 
    C = zeros(size(vi));
    
    sel1 = vs >= vi & vs < vi+1;
    C(sel1) = (1/6)*(vs(sel1)-vi(sel1)).^3;
    
    sel2 = vs >= vi+1 & vs < vi+2;
    C(sel2) = (1/6)*(-3*(vs(sel2)-vi(sel2)-1).^3 + 3*(vs(sel2)-vi(sel2)-1).^2 + 3*(vs(sel2)-vi(sel2)-1)+1);
    
    sel3 = vs >= vi+2 & vs < vi+3;
    C(sel3) = (1/6)*(3*(vs(sel3)-vi(sel3)-2).^3 - 6*(vs(sel3)-vi(sel3)-2).^2 + 4);
    
    sel4 = vs >= vi+3 & vs < vi+4;
    C(sel4) = (1/6)*(1-(vs(sel4)-vi(sel4)-3)).^3;

end