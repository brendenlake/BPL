%
% Given the coarse program structure, enumerate all
% of the attachment relations.
%
% Output:
%  all_R: [ns x 1 cell] each cell contains a cell array
%    of all the relations that stroke could participate in
%
function all_R = cache_enum_all_relations(lib,M)
    all_R = cell(M.ns,1);
    for s=1:M.ns
       prev = 1:s-1; 
       [all_R{s}.list_R,all_R{s}.list_score_type] = enum_all_relations(M.S(prev),lib);
    end
end