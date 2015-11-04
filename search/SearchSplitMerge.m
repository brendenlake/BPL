%
% Greedy search for split/merging strokes
% into sub-strokes
%
function M = SearchSplitMerge(M,lib,verbose)

    assert(M.has_relations());
    
    tabu_item = [];
    score_M = scoreMP(M,lib);
    
    accept = true;
    if verbose, fprintf(1,'  start split/merge search.\n'); end
    while accept
        accept = false;
        
        % get all moves
        [moves,tabu] = get_all_moves(M);
        n = length(moves);
        
        % try all moves
        for i=1:n
            
            if ~isequal(moves{i},tabu_item) % if not tabu item
                Q = do_move(M,moves{i},lib);
                
                score_Q = scoreMP(Q,lib);
                
                if score_Q > score_M % accept or reject
                   M = Q;
                   score_M = score_Q;
                   tabu_item = tabu{i};
                   accept = true;
                   if verbose, print_move(moves{i}); end
                   break
                end
                
            end
        end
        
    end
    if verbose, fprintf(1,'  end split/merge search.\n'); end
    
end

function [moves,tabu] = get_all_moves(M)
    % make a cell array of all possible moves
    [moves1,tabu1] = UtilMP.all_merge_moves(M);
    [moves2,tabu2] = UtilMP.all_split_moves(M);
    moves = [moves1; moves2];
    tabu = [tabu1; tabu2];
end

function Q = do_move(M,move,lib)
    if strcmp(move.type,'split')
        Q = opt_propose_split(M,move.sid,move.bid_start_of_second,lib);
    elseif strcmp(move.type,'merge');
        Q = opt_propose_merge(M,move.i1,lib);
    else
        error('invalid move type'); 
    end
end

function print_move(m)
    % display what move we just did
    if strcmp(m.type,'split')
        fprintf(1,'    split stroke %d at %d\n',m.sid,m.bid_start_of_second);
    else
        fprintf(1,'    merged strokes %d at %d\n',m.i1,m.i2);
    end
end


function Q = opt_propose_split(M,sid,bid_start_of_second,lib)
%
% Propose spliting a stroke in to two sub-components.
%
% Input
%  sid: stroke id
%  bid_start_of_second: sub-stroke id to start new stroke with
%  

    % perform the splitting
    Z = M.S{sid}.copy();
    [S1,S2] = UtilMP.split_stroke(Z,bid_start_of_second);
    Q = M.copy();    
    Q.S =[Q.S(1:sid-1); {S1}; {S2}; Q.S(sid+1:end)];
    
    % clear and optimize relations
    Q.clear_relations();    
    argmax_relations(lib,Q);
    
    % optimize the sub-strokes
    optimize_this_subid(Q,sid,lib);
    optimize_this_subid(Q,sid+1,lib);
end


function Q = opt_propose_merge(M,sid,lib)
%
% Propose merging stroke sid and sid+1
%
% Input
%  sid: stroke id
%  bid_start_of_second: sub-stroke id to start new stroke with
%
  
    % perform the merging
    S1 = M.S{sid}.copy();
    S2 = M.S{sid+1}.copy();
    Q = M.copy();
    Z = UtilMP.merge_strokes(S1,S2);
    Q.S{sid} = Z;
    Q.S(sid+1) = [];
    
    % clear and optimize relations
    Q.clear_relations();    
    argmax_relations(lib,Q);    
    
    % optimize the sub-strokes, gradient, then sub-strokes
    optimize_this_subid(Q,sid,lib);
    argmax_fit_type(Q,lib,sid);
    optimize_this_subid(Q,sid,lib);
end