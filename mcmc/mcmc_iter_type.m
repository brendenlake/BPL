%
% Perform an iteration of MCMC
% for re-sampling all of the type-level
% variables, excluding relation types
%
function mcmc_iter_type(MH,M,lib)

    % shape types
    for sid=1:M.ns
        for bid=1:M.S{sid}.nsub
            MH.mh_shape_type(sid,bid,M,lib);
        end
    end       

    % scale types
    for sid=1:M.ns
        for bid=1:M.S{sid}.nsub
            MH.mh_scale_type(sid,bid,M,lib);
        end
    end

    % specific relation parameters
    for sid=1:M.ns 
       if strcmp(M.S{sid}.R.type,'unihist')
          MH.mh_gobal_position(sid,M,lib);
       elseif strcmp(M.S{sid}.R.type,'mid')
          MH.mh_eval_spot_type(sid,M,lib);
       end
    end

    % sub-stroke ids
    for sid=1:M.ns
        for bid=1:M.S{sid}.nsub
             MH.gibbs_substroke_id(sid,bid,M,lib); 
        end
    end

end