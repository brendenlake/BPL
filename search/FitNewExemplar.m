function [Mfit,score_fit] = FitNewExemplar(I,samples,lib,auto_affine,fit_affine_only)
% FitNewExemplar Fit samples of MotorProgram types 
% to a new image "I" of a new character
%
% Input
%   I : image
%   samples [K x 1 cell]: multiple MotorPrograms as output of MCMC on type-leve
%   parametrs
%   lib : 
%   auto_affine: initialize the affine warp automatically
%   fit_affine_only: fit only the affine warp parameters
%
% Output
%   Mfit: returns an object of type MotorProgramFit   
%
    if ~exist('auto_affine','var')
       auto_affine = true; 
    end
    if ~exist('fit_affine_only','var')
       fit_affine_only = false; 
    end

    Mfit = MotorProgramFit(samples);
    Mfit.I = I;
    
    % initialize the affine warp
    if auto_affine
        UtilMP.set_affine_to_image(Mfit,I);
    end
    
    % blur out the image before starting
    Mfit.blur_sigma = Mfit.parameters.max_blur_sigma;
    Mfit.epsilon = 1e-2;
    
    % do the optimization
    if fit_affine_only % AFFINE
        argmax_fit_affine(Mfit,lib);
        score_fit = scoreMP_fit(Mfit,lib);
    else % REGULAR
        score_fit = argmax_fit_token(Mfit,lib);        
    end
    
end