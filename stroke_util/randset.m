% Pick an index from the vector v,
% where each element is 0 or 1.
% Only pick from the elements denoted by 1.
%
% Output
%  item: scalar indx
function item = randset(v)
    assert(isvector(v));
    assert(all(v==0|v==1));
    indx = find(v);
    z = randint(1,1,[1 length(indx)]);
    item = indx(z);    
end