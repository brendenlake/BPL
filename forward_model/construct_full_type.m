% Construct a full motor program based just on the 
% type-level variables
function Q = construct_full_type(M)
    
    % fill in the structural parameters
    Q = M.copy();
    check_type_complete(Q);
    clear_token_level(Q);
    fill_token_with_type(Q);
    
    % fill in the rendering paramters
    Q.A = [];
    Q.blur_sigma = Q.parameters.min_blur_sigma;
    Q.epsilon = Q.parameters.min_epsilon;
    Q.I = Q.pimg > 0.5;

end

%% make sure that we have type-levels shape parameters
function check_type_complete(Q)
    ns = Q.ns;
    for sid=1:ns
       S = Q.S{sid};
       
       if isempty(S.shapes_type)
          Q.S{sid}.shapes_type = Q.S{sid}.shapes_token;
       end
       
       % check attach relation
       type = Q.S{sid}.R.type;
       if strcmp(type,'mid') && isempty(Q.S{sid}.R.eval_spot_type)
          Q.S{sid}.R.eval_spot_type = Q.S{sid}.R.eval_spot_token; 
       end
       
    end
end

%% make sure that we have nothing at the token level
function clear_token_level(Q)
    ns = Q.ns;
    for sid=1:ns
       Q.S{sid}.pos_token = [];
       Q.S{sid}.shapes_token = [];
       Q.S{sid}.invscales_token = [];
       type = Q.S{sid}.R.type;
       if strcmp(type,'mid')
          Q.S{sid}.R.eval_spot_token = []; 
       end
    end
end

%% fill in the token-level variables with their expectation
function fill_token_with_type(Q)
    ns = Q.ns;
    for sid=1:ns
       
       % fill in token-level relatoins
       type = Q.S{sid}.R.type;
       if strcmp(type,'mid')
          Q.S{sid}.R.eval_spot_token = Q.S{sid}.R.eval_spot_type; 
       end
        
       % fill in the position
       previous_strokes = Q.S(1:sid-1);
       R = Q.S{sid}.R;
       base = getAttachPoint(R,previous_strokes); 
       Q.S{sid}.pos_token = base;
       
       % fill in the shape token
       Q.S{sid}.shapes_token = Q.S{sid}.shapes_type;
       
       % fill in the scale token
       Q.S{sid}.invscales_token = Q.S{sid}.invscales_type;
       
    end
end