%
% Extracts Affine model parameters into a vector.
% Used for optimizing a parse, with fixed type-level parameters,
% to a new image
%
% Extract into a vector:
%  image blur, pixel noise, and affine parameters
% 
function vec_to_model_fit_affine(theta,M)
    
    count = 1;
    
    % pixel noise
    M.epsilon = theta(count);
    count = count + 1;
    
    % image blur
    M.blur_sigma = theta(count);
    count = count + 1;

    % affine warp
    x = count:count+3;
    M.A = vec(theta(x));
    count = count + 4;
  
end