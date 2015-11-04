%
% Resize and center a motor program.
%
% Input
%   M : motor program
%   new_scale : number of pixels spanned by the longest dimension
%     value of < 0 means do not rescale
%  
function my_center_mp(M,new_scale)

    if ~exist('new_scale','var')
        new_scale = 65;       
    end
    
    assert(~isempty(M.I));             
    traj = flatten_nested(M.motor_warped);
    
    curr_rg = range_char(traj);
    if new_scale < 1
        target_rg = curr_rg;
    else
        scale = new_scale / max(curr_rg);   
        target_rg = scale .* curr_rg;
    end        
    imsize = size(M.pimg);                        
    
    UtilMP.set_affine_motor(M, [imsize(1)/2 -imsize(2)/2], target_rg);
    M.I = M.pimg > 0.5;
        
end