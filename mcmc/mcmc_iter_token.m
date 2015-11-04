%
% Perform an iteration of MCMC
% for re-sampling all of the token-level
% variables
%
function mcmc_iter_token(MH,M,lib)

    % shape tokens
    for sid=1:M.ns
        for bid=1:M.S{sid}.nsub
            MH.mh_shape_token(sid,bid,M,lib);
        end
    end       

    % scale tokens
    for sid=1:M.ns
        for bid=1:M.S{sid}.nsub
            MH.mh_scale_token(sid,bid,M,lib);
        end
    end

    % position tokens
    for sid=1:M.ns
       MH.mh_token_position(sid,M,lib); 
    end
    
    % evaluation token spots
    for sid=1:M.ns
       if strcmp(M.S{sid}.R.type,'mid')
          MH.mh_eval_spot_token(sid,M,lib); 
       end        
    end

end