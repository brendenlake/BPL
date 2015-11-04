%
% Replaces model parameters (from a vector) into a the class
% Used for optimizing a parse to an image
%
% Replaces variables:
%  image blur and pixel noise
%  stroke position, scale, and shape
%
% Requirements:
%  shapes_type should be empty, so we can marginalize over it 
%  
function vec_to_model_fit_type(theta,M,list_sid)

    if ~exist('list_sid','var')
       list_sid = 1:M.ns; 
    end

    % Criteria for use
    for i=list_sid
       assert(isempty(M.S{i}.shapes_type));
    end
    
    count = 1;
    
    % pixel noise
    M.epsilon = theta(count);
    count = count + 1;
    
    % image blur
    M.blur_sigma = theta(count);
    count = count + 1;

    for i=list_sid
       
        % stroke position
        x = count:count+1;
        M.S{i}.pos_token = vec(theta(x))';
        count = count + 2;
        
        % sub-stroke shapes
        nn = numel(M.S{i}.shapes_token);
        x = count:count+nn-1;
        M.S{i}.shapes_token = reshape(theta(x),size(M.S{i}.shapes_token));
        count = count + nn;
        
        % sub-stroke scales (type)
        nn = numel(M.S{i}.invscales_type);
        x = count:count+nn-1;
        M.S{i}.invscales_type = reshape(theta(x),size(M.S{i}.invscales_type));
        count = count + nn;
        
        % sub-stroke scales (token)
        nn = numel(M.S{i}.invscales_token);
        x = count:count+nn-1;
        M.S{i}.invscales_token = reshape(theta(x),size(M.S{i}.invscales_token));
        count = count + nn;
        
    end

end