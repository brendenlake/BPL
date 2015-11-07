function [scoreF,score0] = argmax_fit_affine(Mfit_init,lib)
%ARGMAX_FIT_TOKEN Fit a parse to an image, with token-level parameters
%

    % set up objective function
    Mfit = Mfit_init.copy();    
    fmin = @(theta) myscore(theta,Mfit,lib);

    % Set-up the optimization problem  
    [theta0,lb,ub] = model_to_vec_fit_affine(Mfit);    
    options = optimset('Display','off');
    % options = optimset(options,'Display','iter');
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
    vec_to_model_fit_affine(theta,M);
end

% fill the MotorProgram with parameters
% and then score
function minscore = myscore(theta,Mfit,lib)
    Qfit = Mfit.copy(); % we don't want to modify the shared MotorProgramFit base
    refill(theta,Qfit);
    
    PM = Qfit.parameters;
    % score blur, pixel noise affine, and image
    ll = CPD.score_image_blur(Qfit.blur_sigma,PM) + CPD.score_image_noise(Qfit.epsilon,PM);
    ll = ll + CPD.score_affine(lib,Qfit.A) + CPD.score_image(Qfit.I,Qfit.pimg);
    minscore = -ll;
end