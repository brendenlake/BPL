%
% Extracts model parameters into a vector.
% Used for optimizing a parse to an image
%
% Extract into a vector:
%  image blur and pixel noise
%  stroke position, scale, and shape
%
% Requirements:
%  shapes_type should be empty, so we can marginalize over it 
%
% Input
%   list_sid: (default=all) list of stroke-ids which we want to optimize
%    while the others are frozen
function [theta,lb,ub] = model_to_vec_fit_type(M,list_sid)

    if ~exist('list_sid','var')
       list_sid = 1:M.ns; 
    end

    % Criteria for use
    for i=list_sid
       assert(isempty(M.S{i}.shapes_type));
    end

    % count the number of parameters
    ntot = 2; % blur and noise
    for i=list_sid
        ntot = ntot + 2; % for position
        ntot = ntot + numel(M.S{i}.shapes_token);
        ntot = ntot + numel(M.S{i}.invscales_type);
        ntot = ntot + numel(M.S{i}.invscales_token);
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
        
        % sub-stroke scales (type)
        nn = numel(stk.invscales_type);
        x = count:count+nn-1;
        theta(x) = stk.invscales_type;
        lb(x) = 0+ee;
        ub(x) = 1-ee;
        count = count + nn; 
        
        % sub-stroke scales (token)
        nn = numel(stk.invscales_token);
        x = count:count+nn-1;
        theta(x) = stk.invscales_token;
        lb(x) = 0+ee;
        ub(x) = 1-ee;
        count = count + nn;
        
    end

    assert(numel(theta)==ntot);
end