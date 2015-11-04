%
% Center and re-scale all the sub-strokes in a set
% 
% input
%  drawings: [nested cell] each base cell contains a stroke [m x2]
%      it can have an arbitrary nesting structure
%   newscale : length of longest dimension for each stroke
%   
% output
%    drawings_norm : [nested cell] strokes after normalization
%    drawings_centers : [nested cell] previous center of mass for each
%       stroke
%    drawings_scales : [nested cell] previous scale for each stroke
%
%
function [drawings_norm,drawings_centers,drawings_scales] = normalize_dataset(drawings,newscale,verbose)

    if ~exist('verbose','var')
       verbose = true; 
    end
    if verbose, fprintf(1,'get means\n'); end
    drawings_centers = apply_to_nested(drawings,@(stk)com_stk(stk)); % contains mean of each stroke
    if verbose, fprintf(1,'translate\n'); end
    drawings_zeromean = apply_to_nested(drawings,@(stk)offset_stk(stk,com_stk(stk))); % each stroke is now centered
    fscale = @(stk) repmat( newscale ./ max( 1,max( range(stk,1) ) ), [1 2]);
    if verbose, fprintf(1,'get scales\n'); end
    drawings_scales = apply_to_nested(drawings_zeromean, @(stk)fscale(stk)); % contains scale of each stroke
    if verbose, fprintf(1,'rescale\n'); end
    drawings_norm = apply_to_nested(drawings_zeromean, @(stk)rescale_stk(stk,fscale(stk))); % each stroke is now centered and scaled

end