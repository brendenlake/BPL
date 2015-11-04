% function expand_small_strokes
% ---------
% Expand short strokes so they are of the minimum length.
%
% --------
% input
%   dataset: nested cell arrays, where the bottom level
%       is a cell array of "s" strokes. Each stroke
%       is an array of [m x 2] where m varies
%   minlen: the minimum length
%   bool_det: (optional) true if you want deterministic function
function ndataset = expand_small_strokes(dataset,minlen,bool_det)
    if ~exist('bool_det','var')
       bool_det = false; 
    end
    ndataset = apply_to_nested(dataset,@(x)expand(x,minlen,bool_det));
end

function new_stroke = expand(stroke,minlen,bool_det)
    [len,d] = size(stroke);
    assert(d==2); % stroke must have dimensionality two
    if len>= minlen
       new_stroke = stroke;
       return 
    end
    mc = floor(minlen/len);
    assert(mc>=1);
    remainder = minlen - mc*len;
    indx = 1:len;
    newindx = repmat(indx,[1 mc])';
    if bool_det % deterministic
        rperm = 1:len; 
    else
        rperm = randperm(len); 
    end
    remindx = rperm(1:remainder);
    newindx = [newindx; remindx'];
    newindx = sort(newindx);
    assert(length(newindx)==minlen);
    new_stroke = stroke(newindx,:);
end