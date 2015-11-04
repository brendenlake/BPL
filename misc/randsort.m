%
% Just like the "sort" function, 
% except that ties are broken randomly
% rather than preserving the original order
%
% Output
%   B: sorted vector v
%   IX: B = v(IX)
function [B,IX] = randsort(v,direction)

    assert(isvector(v));
    assert(strcmp(direction,'ascend') || strcmp(direction,'descend'));
    
    n = numel(v);
    perm = randperm(n);
    v = v(perm);    
    [B,IX2] = sort(v(:),1,direction);
    
    IX = vec(1:n);
    IX = IX(perm);    
    IX = IX(IX2);    
end