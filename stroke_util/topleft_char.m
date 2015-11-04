%
% Get the top-left point for a given character.
% Only works for characters in MOTOR space 
%
function topleft = topleft_char(char)

    ns = length(char);
    mymax = zeros(ns,2);
    mymin = zeros(ns,2);
    
    % for each stroke
    for i=1:ns
       stk = char{i};
       mymax(i,:) = max_stk(stk);
       mymin(i,:) = min_stk(stk);
    end
    
    % Get overall max and min
    mx = max(mymax,[],1);
    mn = min(mymin,[],1);
    
    % assert(mx(2) <= 0);
     
    % Get the top-left point
    left = mn(1);
    top = mx(2);
    topleft = [left top];    
end