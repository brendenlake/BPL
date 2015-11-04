function [scoreF,score0] = argmax_fit_token(Mfit_init,lib)
% ARGMAX_FIT_TOKEN Fit a parse to an image, with token-level parameters
% only
% 
%  Mfit_init : instance of class MotorProgramFit
%
% Output:
%  scoreF : final score
%  score0 : initial score
% 

    % set up objective function
    Mfit = Mfit_init.copy();    
    fmin = @(theta) myscore(theta,Mfit,lib);

    % Set-up the optimization problem  
    [theta0,lb,ub] = model_to_vec_fit_token(Mfit);    
    options = optimset('Display','off');
    options = optimset(options,'Display','iter');
    options = optimset(options,'TolFun',1e-4,'Algorithm','active-set');

    % Run the optimization
    score0 = -fmin(theta0);
    try
        thetaF = fmincon(fmin,theta0,[],[],[],[],lb,ub,[],options);
    catch
        fprintf(1,'*Warning*: optimization failed. Error caught.\n'); 
        thetaF = theta0;
    end
    scoreF = -fmin(thetaF);
    
    refill(thetaF,Mfit_init);
    
end

function refill(theta,M)
    vec_to_model_fit_token(theta,M);
end

% fill the MotorProgram with parameters
% and then score
function minscore = myscore(theta,Mfit,lib)
    Qfit = Mfit.copy(); % we don't want to modify the shared MotorProgramFit base
    refill(theta,Qfit);
    ll = scoreMP_fit(Qfit,lib);    
    minscore = -ll;
end