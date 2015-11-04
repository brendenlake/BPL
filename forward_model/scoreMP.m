%
% Score a MotorProgram
%   log-probability score as in Equation 1 in the main text
%
% Input
%  M: motor program
%  libclass: object of type LibraryClass
%  
%  By default, it includes all components of the log-score.
%  To choose only specific components, you can use this syntax:
%    'token',true,'type',true,'image',true,'stat',true
%     
%     For the token, type, and image model respectively set to true or
%       false. 'stat' is only needed for alphabet model, which is not 
%       included here and should not be used.
% 
%  Output
%   score : score
%   out : struct that breaks down the score by component
% 
function [score,out] = scoreMP(M,libclass,varargin)

    % should we include type and/or token-level variables?
    assert(isa(libclass,'Library'));
    n = length(varargin);
    if mod(n,2)~=0
       error('the number of optional arguments must be even'); 
    end
    
    defaultinc = (n == 0 || (n==2 && strcmp(varargin{1},'strokes')));
    inc_stat = false;
    if defaultinc
        inc_image = true;
        inc_token = true;
        inc_type = true;
    else % specified in the arguments
        inc_image = nan;
        inc_token = nan;
        inc_type = nan;
    end
    list_sid = 1:M.ns;
    for i=1:2:n-1
       switch varargin{i}
           case 'strokes'
               list_sid = varargin{i+1};
               assert(isnumeric(list_sid));
           case 'token'
               inc_token = varargin{i+1};
               assert(islogical(inc_token));
           case 'type'
               inc_type = varargin{i+1};
               assert(islogical(inc_type));
           case 'image'
               inc_image = varargin{i+1};
               assert(islogical(inc_image));
           case 'stat'
               inc_stat = varargin{i+1};
               assert(islogical(inc_stat));
           otherwise 
               error('invalid input parameters to scoreMP');
       end           
    end
    
    % check that we have specified all the terms to include/exclude
    if ~defaultinc
        if any(isnan([inc_image inc_token inc_type]))
           error('score terms not fully specified');
        end
    end
    PM = M.parameters;    
    out = struct;
    
    % check that we have relations instantiated
    if inc_token || inc_type
        assert(M.has_relations(list_sid));
    end
    
    % type-level score for number of strokes
    if inc_type
        out.number = CPD.score_number(libclass,M.ns);
    end
    
    out.S = cell(M.ns,1);
    for i=list_sid
       MS = M.S{i};
       out.S{i} = struct;
       
       % type-level scores
       if inc_type
           out.S{i}.sequence = CPD.score_sequence(libclass,M.ns,MS.ids);          
           out.S{i}.scales_type = CPD.score_invscale_type(libclass,MS.invscales_type,MS.ids);
           out.S{i}.relation = CPD.score_relation_type(libclass,MS.R);
           if ~isempty(MS.shapes_type)
              out.S{i}.shape_type = CPD.score_shape_type(libclass,MS.shapes_type,MS.ids); 
           end
       end
       
       % token-level scores
       if inc_token
           if isempty(MS.shapes_type)
               out.S{i}.shape_marginalize = CPD.score_shape_marginalize(libclass,MS.shapes_token,MS.ids);
           else
               out.S{i}.shape_token = CPD.score_shape_token(libclass,MS.shapes_token,MS.shapes_type);
           end
           out.S{i}.scales_token = CPD.score_invscale_token(libclass,MS.invscales_token,MS.invscales_type);
           out.S{i}.position = CPD.score_position(libclass,MS.pos_token,MS.R,M.S(1:i-1));
           if strcmp(MS.R.type,'mid')
               if isempty(MS.R.eval_spot_type)
                   out.S{i}.relation_token = CPD.score_relation_token_approx_marginalize(libclass,MS.R.eval_spot_token);                   
               else
                   out.S{i}.relation_token = CPD.score_relation_token(libclass,MS.R.eval_spot_token,MS.R.eval_spot_type);
               end
           end
       end       
       
    end
    
    % token-level image variables
    if inc_token
        out.image_blur = CPD.score_image_blur(M.blur_sigma,PM);
        out.image_noise = CPD.score_image_noise(M.epsilon,PM);
        out.A = CPD.score_affine(libclass,M.A);
    end
    
    % include the stat penalty term
    if inc_stat
        error('This term should not be used in the score.');
        out.stat = CPD.score_stat(libclass,M);
    end
    
    % score the image
    if inc_image
        out.image = CPD.score_image(M.I,M.pimg);
    end
    
    score = sumfields(out);
    if isnan(score)
       assert(false); 
    end
end

function total = sumfields(S)
%
% Add up all the scalar, numeric fields
% in a structure S. Recursively call
% on cell arrays of fields
%
    if isempty(S)
       total = 0;
       return
    end
    assert(isstruct(S));
    total = 0;
    fds = fieldnames(S);
    for i=1:length(fds)
       field = S.(fds{i});
       if ~isempty(field)
           
          if isnumeric(field) % if number
              total = total + sum(field(:));
          
          elseif iscell(field) % if cell, call recursively
              for j=1:length(field)
                total = total + sumfields(field{j});
              end
              
          else
              error('invalid field type'); 
          end
          
       end
    end
end