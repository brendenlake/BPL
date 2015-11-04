function [scoreF,score0] = argmax_fit_type(Minit,lib,list_sid,verbose)
% ARGMAX_FIT_TYPE... Fit a parse to an image, with type-level parameters
% only
%
% Input
%   Minit : instance of MotorProgram class
%   list_sid: (default=all) list of stroke-ids which we want to optimize
%    while the others are frozen
%   verbose: (true/false) describe algorithm progress
%
% Output:
%  scoreF : final score
%  score0 : initial score
% 

    M = Minit.copy();
    if ~exist('list_sid','var')
       list_sid = 1:M.ns; 
    end
    if ~exist('verbose','var')
       verbose = false;
    end
    
    has_relations = M.has_relations(list_sid);
    if has_relations
        fmin = @(theta) myscore_HasRel(theta,M,lib,list_sid);
    else
        all_R = cache_enum_all_relations(lib,M);
        fmin = @(theta) myscore_NoRel(theta,M,lib,list_sid,all_R);
    end

    % Set-up the optimization problem  
    [theta0,lb,ub] = model_to_vec_fit_type(M,list_sid);    
    options = optimset('Display','off');
    if verbose
        options = optimset(options,'Display','iter');
    end
    options = optimset(options,'TolFun',1e-4,'Algorithm','active-set');
     
    % Run the optimization
    %tic
    score0 = -fmin(theta0);   
    try
        thetaF = fmincon(fmin,theta0,[],[],[],[],lb,ub,[],options);
    catch
        fprintf(1,'*Warning*: optimization failed. Error caught.\n'); 
        thetaF = theta0;
    end
    scoreF = -fmin(thetaF);
    %t = toc;
    
    refill(thetaF,Minit,list_sid);
    
    % make sure we have not added/removed relations
    assert(Minit.has_relations==has_relations);
    
end

function refill(theta,M,list_sid)
    vec_to_model_fit_type(theta,M,list_sid);
end

% fill the MotorProgram with parameters
% and then score
function minscore = myscore_HasRel(theta,M,lib,list_sid)
    Q = M.copy(); % we don't want to modify the shared MotoProgram base
    refill(theta,Q,list_sid);
    ll = scoreMP(Q,lib,'strokes',list_sid,'type',true,'token',true,'stat',false,'image',true);    
    minscore = -ll;
end

% fill the MotorProgram with parameters
% and then optimize the relations.
% Finally, return the score "minscore"
function minscore = myscore_NoRel(theta,M,lib,list_sid,all_R)
    Q = M.copy(); % we don't want to modify the shared MotoProgram base
    refill(theta,Q,list_sid);
    argmax_relations(lib,Q,all_R,list_sid);
    ll = scoreMP(Q,lib,'strokes',list_sid,'type',true,'token',true,'stat',false,'image',true); 
    minscore = -ll;
end