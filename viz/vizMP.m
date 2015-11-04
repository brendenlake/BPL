%
% Visualize a MotorProgram
%
% additional arguments 'motor' plots the motor trace,
% and the arguments 'image' and 'pmap' can select
% betwen the image and the probability map
%
%  Other parameters
%  'size_start' : where 0 if we don't want to plot the start points
%  'lw' : line width for strokes
function vizMP(M,varargin)

    assert(isa(M,'MotorProgram'));
   
    n = length(varargin);
    use_pimg = any(strcmp('pmap',varargin)) || any(strcmp('pimg',varargin));
    use_img = any(strcmp('image',varargin));
    assert(~(use_pimg && use_img));
    use_motor = any(strcmp('motor',varargin));
    
    % what is the size of the start marker?
    fstart = find(strcmp('size_start',varargin));
    size_start = [];
    if ~isempty(fstart)
        size_start = varargin{fstart+1};
    end
    assert(isnumeric(size_start));
    
    % what is the line width?
    flw = find(strcmp('lw',varargin));
    lw = [];
    if ~isempty(flw)
       lw = varargin{flw+1}; 
    end

    
    if use_motor
        if use_pimg            
            plot_motor_to_image(M.pimg,M.motor_warped,size_start,lw);
        else                       
            plot_motor_to_image(M.I,M.motor_warped,size_start,lw);            
        end
    else
        if use_pimg
            plot_image_only(M.pimg); 
        else
            plot_image_only(M.I); 
        end
    end
            
end