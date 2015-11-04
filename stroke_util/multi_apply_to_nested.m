function Ac = multi_apply_to_nested(fnc,nested1,nested2,nested3)
    % apply function handle fnc(item1,item2,item3) where items
    % are the base items in each collection of nested elements

    nopt = nargin-1;
    if nopt == 1
        Ac = process(fnc,nested1);
    elseif nopt == 2
        Ac = process(fnc,nested1,nested2);
    elseif nopt == 3;
        Ac = process(fnc,nested1,nested2,nested3);
    else
        error('invalid number of arguments');
    end
   
end

% Recursive function 
function Ac = process(fnc,A,B,C)

    % Assign varialbes
    nopt = nargin - 1;
    if nopt < 2
        B = [];
    end
    if nopt < 3
        C = [];
    end
    
    isbase = @(x) ~iscell(x) && ~isempty(x);
    
    % Base case
    if (isbase(A) || isbase(B) || isbase(C))
       %assert(~iscell(B) && ~iscell(C))
       if nopt == 1
           Ac = fnc(A);
       elseif nopt == 2
           Ac = fnc(A,B);
       elseif nopt == 3
            Ac = fnc(A,B,C);
       else
           error('invalid number of arguments');
       end
       
    % Recursive case
    else
       n = numel(A);
       Ac = A;
       for i=1:n
            if nopt == 1
                Ac{i} = process(fnc,A{i});
            elseif nopt == 2
                Ac{i} = process(fnc,A{i},B{i});
            elseif nopt == 3
                Ac{i} = process(fnc,A{i},B{i},C{i});
            else
                error('invalid number of arguments');
            end
       end
    end
end