%
% Extracts Affine model parameters into a vector.
% Used for optimizing a parse, with fixed type-level parameters,
% to a new image
%
% Extract into a vector:
%  image blur, pixel noise, and affine parameters
%
function [theta,lb,ub] = model_to_vec_fit_affine(M)

    % count the number of parameters
    ntot = 2; % blur and noise
    ntot = ntot + 4; % affine parameters

    theta = zeros(ntot,1);
    lb = -inf(ntot,1);
    ub = inf(ntot,1);
    count = 1;
    
    % pixel noise    
    theta(count) = M.epsilon;
    lb(count) = M.parameters.min_epsilon;
    ub(count) = M.parameters.max_epsilon;
    count = count + 1;
    
    % image blur
    theta(count) = M.blur_sigma;
    lb(count) = M.parameters.min_blur_sigma;
    ub(count) = M.parameters.max_blur_sigma;
    count = count + 1;
    
    % affine warp
    x = count:count+3;
    theta(x) = M.A;
    max_scale = M.parameters.max_affine_scale_change;
    max_shift = M.parameters.max_affine_shift_change;
    lb(x) = [1/max_scale; 1/max_scale; -max_shift; -max_shift];
    ub(x) = [max_scale; max_scale; max_shift; max_shift];
    count = count + 4;

    assert(numel(theta)==ntot);
end