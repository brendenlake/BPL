function optimize_this_subid(Q,sid,lib,verbose)
    % optimize all of the sub-strokes in stroke number sid
    %
    % Input
    %  Q : MotorProgram
    %  sid : stroke id
    % 
    if ~exist('verbose','var')
       verbose = false; 
    end
    
    
    assert(numel(sid)==1);
    iter = 1;
    nsub = Q.S{sid}.nsub;
    bool_change = true(nsub,1);
    if verbose, fprintf(1,'('); end
    while any(bool_change) && iter < 10
        for bid=1:nsub % each substroke
            score = score_all_subid(Q,sid,bid,lib);
            curr_id = Q.S{sid}.ids(bid);
            new_id = argmax(score);
            Q.S{sid}.ids(bid) = new_id;
            bool_change(bid) = curr_id ~= new_id;
        end       
        if nsub==1, bool_change = false; end % no need to repeat
        if verbose, fprintf(1,'%d, ',iter); end
        iter = iter + 1;
    end
    if verbose, fprintf(1,'iter.)\n'); end
end