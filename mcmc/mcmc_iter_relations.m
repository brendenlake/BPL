%
% Perform an iteration of MCMC
% for re-sampling the relations, which 
% blends type-level and token-levle variables
%
function mcmc_iter_relations(MH,M,lib)

    % relations
    for sid=2:M.ns
        MH.mh_relation(sid,M,lib);
    end
    
    % evaluation token spots
    for sid=1:M.ns
       if strcmp(M.S{sid}.R.type,'mid')
          MH.mh_eval_spot_token(sid,M,lib); 
       end        
    end

end