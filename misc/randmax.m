%
% Find the maximum element of a vector v,
% where ties are broken randomly.
%
% Output
%   C: [scalar] maximum value
%   I: [scalar] index in v, where ties are broken randomly
%
function [C,I] = randmax(v)
    assert(isvector(v));    
    C = max(v);
    maxidx = find(v==C);
    I = rand_discrete(maxidx);
end