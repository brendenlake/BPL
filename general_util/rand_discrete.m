%
% Sample from a discrete distribution on the
% elements of vcell
%
% Input
%  vcell: [n x 1 cell] or [n x 1] vector
%  wts: [n x 1] weights (unnormzlied), default = uniform
%
% Output
%  item: selected item
%  indx: index in original set
%
function [item,indx] = rand_discrete(vcell,wts)

    assert(isvector(vcell));
    n = numel(vcell);

    if ~exist('wts','var')
       wts = ones(n,1) ./ n;        
    else
       assert(numel(wts)==n); 
       wts = wts ./ sum(wts); 
    end

    indx = find(mnrnd(1,wts));
    if numel(indx)>1
       wts 
       error('invalid sample from discrete distribution'); 
    end
    
    if iscell(vcell)
        item = vcell{indx};
    else
        item = vcell(indx); 
    end
end