% Flatten a nested-cell array structure
% into a single cell-array
%
% Input
%  dataset: nested cell-array structure
%
% Output
%  C: flattened-structure
%
function C = flatten_nested(dataset)
    
    num_el = count_elements(dataset);
    C = take_elements(dataset,num_el);
end

% Count the number of elements
% in the nested-array structure
function count = count_elements(A)

    count = 0;
    find_bottom(A);
    
    % Recursive helper for counting
    function find_bottom(A)
        if ~iscell(A)
           count = count + 1;
        else
           n = numel(A);
           for i=1:n
             find_bottom(A{i});
           end
        end
    end

end

% Count the number of elements
% in the nested-array structure
function C = take_elements(A,num_el)
    
    C = cell(num_el,1);

    count = 1;
    get_bottom(A);
    
    % Recursive helper for counting
    function get_bottom(A)
        if ~iscell(A)
           C{count} = A;
           count = count + 1;
        else
           n = numel(A);
           for i=1:n
             get_bottom(A{i});
           end
        end
    end

end