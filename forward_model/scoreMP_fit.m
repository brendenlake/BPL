function score = scoreMP_fit(Mfit,lib)
    % 
    % From a set of samples of type-level parameters (Mfit), score
    % the probability of a new image Mfit.I
    %
    % Computes the following score:
    % 
    %   log P(I(T) | \theta(T)) + logsumexp(log P(\theta(T) | \phi[j] )) - log(N)
    %     for a small set j=1,..,N
    %
    %  which is the inner summation of Equation S10
    %
    nsamp = Mfit.nsamples;
    assert(Mfit.isEmptyType); % ensures we don't modify the the 
        % argument while scoring it
    
    % score all token varibles that depend on the type-level parameters
    score_type = zeros(nsamp,1);
    for i=1:nsamp
       MType = Mfit.setSampleType(i); 
       score_type(i) = scoreMP(MType,lib,'type',false,'token',true,'image',false);
    end
    Mfit.clearSampleType();
    
    score_type = score_type - log(nsamp);    
    marginal_score = logsumexp(score_type);
    
    % score the image, which depends only on token-level parameters
    score_image = scoreMP(Mfit,lib,'image',true,'type',false,'token',false);
    
    % compute the combined score
    score = score_image + marginal_score;
    
end