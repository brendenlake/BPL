%
% Replaces model parameters (from a vector) into an object.
% Used for optimizing a parse, with fixed type-level parameters,
% to a new image
%
% Replace variables:
%  image blur, pixel noise, and affine parameters
%  stroke position, token scale, and token shape
%  
function vec_to_model_fit_token(theta,M,list_sid)

    if ~exist('list_sid','var')
       list_sid = 1:M.ns; 
    end
    
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
    
    for i=list_sid
       
        % stroke position
        x = count:count+1;
        M.S{i}.pos_token = vec(theta(x))';
        count = count + 2;
        
        % sub-stroke shapes (token)
        nn = numel(M.S{i}.shapes_token);
        x = count:count+nn-1;
        M.S{i}.shapes_token = reshape(theta(x),size(M.S{i}.shapes_token));
        count = count + nn;
        
        % sub-stroke scales (token)
        nn = numel(M.S{i}.invscales_token);
        x = count:count+nn-1;
        M.S{i}.invscales_token = reshape(theta(x),size(M.S{i}.invscales_token));
        count = count + nn;
        
        % eval_spot_token
        if ~isempty(M.S{i}.R) && strcmp(M.S{i}.R.type,'mid')
            M.S{i}.R.eval_spot_token = theta(count);
            count = count + 1;
        end        
    end

end