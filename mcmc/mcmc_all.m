function samples = mcmc_all(Q,lib,nsamp,varargin)
    % MCMC for computing local parameter variance.
    %
    % Input
    %  nsamp: the number of samples
    %
    % Include these string parameters to resample
    %   type
    %   token
    %   relations
    % 
    %  By default these variables are not resampled 
    %
    % Optional input strings
    %  debug: use full score for debuggin
    %  verbose: track sampling as we go
    %
 
    narg = numel(varargin);
    verbose = false;
    sample_type = false;
    sample_relations = false;
    sample_token = false;
    debug = false;
    for i=1:narg
        switch varargin{i}
            case 'verbose'
                verbose = true;
            case 'debug'
                debug = true;
            case {'relation','R','relations'}
                sample_relations = true;
            case 'token'
                sample_token = true;
            case 'type'
                sample_type = true;
            otherwise
                error('invalid input parameter');
        end
    end
    
    if ~(sample_type || sample_token || sample_relations)
       error('empty sampler'); 
    end

    MH = MCMC(debug);

    M = Q.copy();
    samples = cell(nsamp,1);
    samples_score = nan(nsamp,1);
    
    % for each sample
    for is = 1:nsamp 
        
        if sample_type
           mcmc_iter_type(MH,M,lib);
        end
        
        if sample_token
           mcmc_iter_token(MH,M,lib); 
        end
        
        if sample_relations
           mcmc_iter_relations(MH,M,lib); 
        end
          
        samples{is} = M.copy();
        samples_score(is) = scoreMP(M,lib);
   
        assert(~isinf(samples_score(is)));
        
        % visualization
        if verbose
            fprintf(1,'Sample %d\n',is);
            figure(1);
            clf
            xlabel('iterations');
            ylabel('log-score');
            plot(1:is,samples_score(1:is));
            drawnow
            pause(0.1);
        end
                
    end
    
    if verbose
        MH.acceptance_ratios();
    end

end