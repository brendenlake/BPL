% Apply a operator (fnc) to
% each of the strokes in a character.
%
% Input
%   char: [ns x 1 cell] array of strokes
%
% Output
%   char: [ns x 1 cell] array of strokes
%
function char = apply_each_stroke(char,fnc)
    assert(iscell(char));
    assert(~iscell(char{1}));
    char = apply_to_nested(char,fnc);
end