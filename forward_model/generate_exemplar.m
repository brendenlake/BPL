%
% Stochastic function that generates new examples from the prior,
% given a type-level specification from generate_character().
% 
% Input
%  template: MotorProgram class with R,ids,shapes_type,invscales_types
%    defined 
%
% Output
%  M: fully-specified MotorProgram
function [M,template] = generate_exemplar(template,libclass)

    M = MotorProgram(template);
    
    % sample stroke parameters
    for i=1:M.ns
        if strcmp(M.S{i}.R.type,'mid')
            M.S{i}.R.eval_spot_token = CPD.sample_relation_token(libclass,M.S{i}.R.eval_spot_type);
        end        
        M.S{i}.pos_token = CPD.sample_position(libclass,M.S{i}.R,M.S(1:i-1));        
        M.S{i}.shapes_token = CPD.sample_shape_token(libclass,M.S{i}.shapes_type);
        M.S{i}.invscales_token = CPD.sample_invscale_token(libclass,M.S{i}.invscales_type);
    end
    
    % sample affine warp
    M.A = CPD.sample_affine(libclass);  
    
    % set rendering parameters to minimum noise
    M.blur_sigma = template.parameters.min_blur_sigma;
    M.epsilon = template.parameters.min_epsilon;
    
    % sample rendering parameters
    %M.blur_sigma = CPD.sample_image_blur(template.parameters);
    %M.epsilon = CPD.sample_image_noise(template.parameters);
    
    % sample the image
    M.I = CPD.sample_image(M.pimg);
    
end