%
% Make a string by combining a sequence of
% both strings and numbers (rounded), where each
% number is converted to a string
%
function str = makestr(varargin)

    n = numel(varargin);
    str = [];
    for i=1:n       
        A = varargin{i};
        if ischar(A)
            str = [str A];
        elseif isnumeric(A)
            str = [str num2str(A,3)];
        else
            error('invalid input type'); 
        end        
    end
    
end