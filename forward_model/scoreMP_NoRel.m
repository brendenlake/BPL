%
% Score the MotorProgram,
%   when the stroke relations are undefined.
%   We implicity maximize over relations, but
%   but return an unmodified object.
%
% Input
%  M: motor program
%  lib: library
%  
%  see scoreMP for description of optional parameters
%
function [score,out] = scoreMP_NoRel(M,lib,varargin)
    
    % make sure we don't have relations for the relevant strokes
    list_sid = 1:M.ns;
    for i=1:2:length(varargin)-1
       switch varargin{i}
           case 'strokes'
               list_sid = varargin{i+1};
               assert(isnumeric(list_sid));
       end           
    end
    assert(~M.has_relations(list_sid));
   
    Q = M.copy();
    argmax_relations(lib,Q,[],list_sid);
    [score,out] = scoreMP(Q,lib,varargin{:});
end