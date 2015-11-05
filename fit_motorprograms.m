%
% Find a set of motor programs for an image
%
% Input
%   I: binary image, where true is "black" pixel
%   K: number 
%   verbose: (true/false) report progress and visualize results?
%   include_mcmc: (true/false) estimate local variance at type-level? 
%       (this is needed for exemplar generation and classification)
%   fast_mode: (true/false)
%   
% Output
%   G: structure to store output
function G = fit_motorprograms(I,K,verbose,include_mcmc,fast_mode)

    ps = defaultps;
    if ~exist('K','var') || isempty(K) 
       K = ps.K; 
    end
    if ~exist('verbose','var')
       verbose = false; 
    end
    if ~exist('include_mcmc','var')
       include_mcmc = true; 
    end
    if ~exist('fast_mode','var')
       fast_mode = false;
    end
    
    load(ps.libname,'lib');
    lib.ncpt;
    
    % generate parses
    initMP = generate_random_parses(I,lib,K,verbose);
    nl = numel(initMP);
    init_scores = zeros(nl,1);
    for i=1:nl
        init_scores(i) = scoreMP_NoRel(initMP{i},lib);
    end      
    
    % visualize initial parses
    if verbose
        sz = [500 500]; % figure size      
        h = figure;
        pos = get(h,'Position');
        pos(3:4) = sz;
        set(h,'Position',pos);
        nrow = ceil(sqrt(nl));        
        for i=1:nl
            subplot(nrow,nrow,i);
            vizMP(initMP{i},'motor');
            lb = num2str(init_scores(i),4);
            if i==1, lb = ['initial score: ' lb]; end
            title(lb);
        end
        pause(0.1);
        drawnow
    end

    % run search for each candidate
    finalMP = cell(nl,1);
    for i=1:nl
        fprintf(1,'\nOptimizing parse %d of %d\n',i,nl);
        finalMP{i} = SearchForParse(initMP{i},lib,verbose,fast_mode);
    end
    final_scores = zeros(nl,1);
    for i=1:nl
       final_scores(i) = scoreMP(finalMP{i},lib); 
    end
    
    % sort programs from best to worst
    [~,score_indx] = sort(final_scores,1,'descend');
    finalMP = finalMP(score_indx);
    final_scores = final_scores(score_indx);
    
    % visualize the final parses
    if verbose
        sz = [500 500]; % figure size      
        h = figure;
        pos = get(h,'Position');
        pos(3:4) = sz;
        set(h,'Position',pos);
        nrow = ceil(sqrt(nl));        
        for i=1:nl
            subplot(nrow,nrow,i);
            vizMP(finalMP{i},'motor');
            lb = num2str(final_scores(i),4);
            if i==1, lb = ['final score: ' lb]; end
            title(lb);
        end
        pause(0.1);
        drawnow
    end
        
    % MCMC to estimate local variance at the type level
    samples_type = [];
    if include_mcmc
        fprintf(1,'\nEstimating local type-level variance for parse...\n',i,nl);
        for j=1:nl
            fprintf(1,'%d,',j);
            samplesM = RunMCMCType(finalMP{j},lib);
            nsamp = numel(samplesM);
            for i=1:nsamp
                samplesM{i}.lightweight; % reduce memory size
            end
            samples_type{j} = samplesM;
        end
        fprintf(1,'done.\n');
        samples_type = samples_type(:);
    end
    
    % return structure
    G = struct;
    G.models = finalMP;
    G.scores = final_scores;
    G.samples_type = samples_type;
    G.img = I;
    G.PM = ps;
    
end