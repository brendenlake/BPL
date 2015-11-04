%
% Extracts model parameters into a vector.
% Used for optimizing a parse, with fixed type-level parameters,
% to a new image
%
% Extract into a vector:
%  image blur, pixel noise, and affine parameters
%  stroke position, token scale, and token shape
%
% Requirements:
%  everything should be instantiated 
%
% Input
%   list_sid: (default=all) list of stroke-ids which we want to optimize
%    while the others are frozen
function [theta,lb,ub] = model_to_vec_fit_token(M,list_sid)

    if ~exist('list_sid','var')
       list_sid = 1:M.ns; 
    end

    % count the number of parameters
    ntot = 2; % blur and noise
    ntot = ntot + 4; % affine parameters
    for i=list_sid
        ntot = ntot + 2; % for position
        ntot = ntot + numel(M.S{i}.shapes_token);        
        ntot = ntot + numel(M.S{i}.invscales_token);
        R = M.S{i}.R;
        if ~isempty(R) && strcmp(R.type,'mid')
            ntot = ntot + 1;
        end
    end

    theta = zeros(ntot,1);
    lb = -inf(ntot,1);
    ub = inf(ntot,1);
    count = 1;
    
    % margin of error for inequalities
    ee = 1e-4;
    
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
        
    for i=list_sid
        
        stk = M.S{i};
        
        % stroke position
        x = count:count+1;
        theta(x) = stk.pos_token;
        lb(x) = [0; -M.parameters.imsize(2)]+ee;
        ub(x) = [M.parameters.imsize(1); 0]-ee;
        count = count + 2;
        
        % sub-stroke shapes (token)
        nn = numel(stk.shapes_token);
        x = count:count+nn-1;
        theta(x) = stk.shapes_token;
        count = count + nn;
        
        % sub-stroke scales (token)
        nn = numel(stk.invscales_token);
        x = count:count+nn-1;
        theta(x) = stk.invscales_token;
        lb(x) = 0+ee;
        ub(x) = 1-ee;
        count = count + nn;
        
        % eval_spot_token
        if ~isempty(stk.R) && strcmp(stk.R.type,'mid')
            theta(count) = stk.R.eval_spot_token;
            ncpt = size(stk.shapes_token,1);
            [~,lb(count),ub(count)] = bspline_gen_s(ncpt,1);
            count = count + 1;
        end
    end

    assert(numel(theta)==ntot);
end