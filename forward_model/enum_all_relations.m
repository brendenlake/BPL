function [list_R,list_score_type] = enum_all_relations(previous_strokes,lib)
    % enumerate all possible ways a new stroke
    % can relate to the previous strokes
    %
    % This instantiates the discrete parameters, but not the  
    % continuous ones like continuous position and attachment.
    %
    % Input
    %  previous_strokes: cell array of stroke objects
    %  
    % Output
    %   list_R: cell array of possible relations
    %   list_score_type : [n x 1] cached scores from
    %           CPD.score_relation_type, where "nan" indicates 
    %           a failure to cache (only the "unihist" case)
    ncpt = lib.ncpt;
    nprev = length(previous_strokes);    
    list_R{1} = RelationIndependent('unihist',nprev);
    for prev=1:nprev
       list_R = [list_R; {RelationAttach('start',nprev,prev)}];
       list_R = [list_R; {RelationAttach('end',nprev,prev)}];       
       for b=1:previous_strokes{prev}.nsub
            nsub = previous_strokes{prev}.nsub;
            list_R = [list_R; {RelationAttachAlong('mid',nprev,prev,nsub,b,ncpt)}];
       end
    end
    
    list_score_type = nan(size(list_R));
    for i=2:length(list_R)
        list_score_type(i) = CPD.score_relation_type(lib,list_R{i});
    end
    
end