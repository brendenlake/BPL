function M = SearchForParse(M,lib,verbose,fast_mode)
    % Search algorithm for finding a good parse
    %
    % Input
    %  M : MotorProgram
    %  verbose: display steps as we go?
    %  fast_mode: (true/false) no gradient search (by far the slowest part)

    if ~exist('verbose','var')
       verbose = false;
    end
    if ~exist('fast_mode','var')
       fast_mode = false; 
    end
    
    if isempty(M)
       return 
    end
    
    if verbose, fprintf(1,'searching for parse (BEGINS)\n'); end
    Do = SearchMoves(M,lib,verbose,fast_mode);   
    Do.disp_score();
    
    % gradient search
    Do.move_opt_grad();    
    
    % The stroke sub-ids have likely changed after optimization. If we 
    % don't update them now, we will be giving an advantage to stroke
    % flipping.
    Do.move_opt_subids();
    
    % gradient search
    Do.move_opt_grad();
    Do.disp_score();
    
    % optimize the direction, order, and relations between strokes
    Do.move_opt_dir_order_rel();
    
    % run a split/merge search
    Do.move_split_merge();
    
    M = Do.M; 
    if verbose, fprintf(1,'searching for parse (ENDS)\n'); end
end