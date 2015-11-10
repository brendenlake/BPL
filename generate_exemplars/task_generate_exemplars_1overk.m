%
% Generate new exemplars.
%
% Rather than re-sampling from the original log-probability scores,
% take the rank score and sample according to 1/sigma(i)
% where sigma(i) is the rank score from highest to lowest score.
%
% Input
%  G: result of model fitting
%  lib: library
%  nsamp: (default=1) number of samples
%
% Output
%  samples: [nsamp x 1 cell] samples
%  types: [nsamp x 1 cell] samples before token-level is resampled
%
function [samples,types] = task_generate_exemplars_1overk(G,lib,nsamp)

    if ~exist('nsamp','var')
       nsamp = 1; 
    end
    
    % re-score the hypotheses by their rank score instead
    % of the real score
    wts = rescore_by_rank(G.scores);
    G.scores = log(wts);
    
    [samples,types] = task_generate_exemplars(G,lib,nsamp);
end

%
% Input
%  G: mixture model struct
%  lib: library
%  nsamp: (default=1) number of samples
%
% Output
%  samples: [nsamp x 1 cell] samples
%  types: [nsamp x 1 cell] samples before token-level is resampled
%
function [samples,types] = task_generate_exemplars(G,lib,nsamp)

    if ~exist('nsamp','var')
       nsamp = 1; 
    end

    % probability of choosing each parse
    logwts = G.scores;
    wts = exp(logwts-logsumexp(logwts(:)));
    samples = cell(nsamp,1);
    types = cell(nsamp,1);
    for iter = 1:nsamp
        
        % choose the parse
        [M,kindx] = rand_discrete(G.models,wts);
        
        % choose the type-level resampling
        Q = rand_discrete(G.samples_type{kindx});
        Q = Q.copy();
        Q.I = M.I;
        
        % resample at the token-level
        types{iter} = Q;
        Q = generate_exemplar(Q.copy(),lib);
        samples{iter} = Q;
    end

end

%
% Transform a score into 1/sigma(i), where
% sigma(i) is the rank (from highest to lowest)
% of the score.
%
% Input
%  scores: [K x 1] vector of scores where higher is better
%
% Output
%  wts: [K x 1] new scores (not in log space) 
%
function wts = rescore_by_rank(scores)

    assert(isvector(scores));
    K = numel(scores);
    
    [~,rank_indx] = sort(scores(:),1,'descend');
    
    wts = zeros(K,1);
    for i=1:K
        wts(i) = 1 ./ rank_indx(i); 
    end
    logwts = log(wts);
    wts = exp(logwts - logsumexp(logwts(:)));
    
    % Assure the precise normalization of the weights 
    sum_wts = sum(wts);
    assert(aeq(sum_wts,1));
    wts = wts ./ sum(wts);  

end