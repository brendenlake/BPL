%
% Classify a trajectory as one of the primitives in the library
%
% Input
%  traj: [n x 2] trajectory
%  scale: [scalar] stroke scale fator
%  
% Output
%
%  structure with fields
%  .indx: best fitting component
%  .bpsline: [ncpt x 2] spline
%  .score: [scalar]
%
function a = classify_traj_as_subid(traj,scale,lib)

    if numel(scale)==2
       scale = scale(1); 
    end
    assert(isscalar(scale));
    [ncat,dim] = size(lib.shape.mu);
    ncpt = dim/2;
    minlen = ncpt*2;

    % fit a spline
    bool_det = true; % determinstic?
    celltraj = expand_small_strokes({traj},minlen,bool_det);
    traj = celltraj{1};
    bspline = fit_bspline_to_traj(traj,ncpt);
    
    rep_shape_type = repmat(bspline,[1 1 ncat]);
    invscale = 1./scale;
    rep_invscale_type = repmat(invscale,[ncat 1]);
    
    % score the various model components
    shape_score = CPD.score_shape_marginalize(lib,rep_shape_type,1:ncat);
    scale_score = CPD.score_invscale_type(lib,rep_invscale_type,1:ncat);
    score = shape_score + scale_score + log(lib.shape.mixprob(:));
    
    % pick the best component
    [best_score,windx] = max(score);
    a = struct;
    a.indx = windx;
    a.bspline = bspline;
    a.score = best_score;
end