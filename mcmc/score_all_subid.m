function score = score_all_subid(M,sid,bid,lib)
%
% For a MotorProgram M, compute the log-score of switching
% its "sub-stroke id" to every other possibility
%
% Input
%  sid: stroke index
%  bid: sub-stroke index
%
% Output
%  score: [nclust x 1] score for each possibility, using all CPDs
%     that directly depend on it
%
    assert(numel(sid)==1);
    assert(numel(bid)==1);
    ncat = lib.N;
    
    seq_score = zeros(ncat,1);
    
    % Score the transition TO the current primitive   
    if bid == 1 % first primitive
        seq_score = seq_score + lib.logStart(:);
    else % not the first
        prev_subid = M.S{sid}.ids(bid-1);
        seq_score = seq_score + vec(lib.logT(prev_subid,1:ncat));
    end
    
    % Score the transition AWAY from the current primitive
    if bid < M.S{sid}.nsub
        next_subid = M.S{sid}.ids(bid+1);
        seq_score = seq_score + vec(lib.logT(1:ncat,next_subid));
    end
    
    % Scale term
    invscale_type = M.S{sid}.invscales_type(bid);
    rep_invscale_type = repmat(invscale_type,[ncat 1]);
    
    % Shape term
    if ~isempty(M.S{sid}.shapes_type)
        shape_type = M.S{sid}.shapes_type(:,:,bid);
        rep_shape_type = repmat(shape_type,[1 1 ncat]);
        shape_score = CPD.score_shape_type(lib,rep_shape_type,1:ncat);
    else
        shape_token = M.S{sid}.shapes_token(:,:,bid);
        rep_shape_token = repmat(shape_token,[1 1 ncat]);
        shape_score = CPD.score_shape_marginalize(lib,rep_shape_token,1:ncat);
    end
    
    % compute scores        
    scale_score = CPD.score_invscale_type(lib,rep_invscale_type,1:ncat);
    score = seq_score + scale_score + shape_score;
end