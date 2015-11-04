%
% Find a set of probable parses for image I
%
% Input
%   I: binary image, where true is "black" pixel
%   K: number 
% 
% Output
%   G: structure to store output
function G = parser(I,K,verbose)

    ps = defaultps;
    if ~exist('K','var') || isempty(K) 
       K = ps.K; 
    end
    if ~exist('verbose','var')
       verbose = false; 
    end
    
    load(ps.libname,'lib');
    lib.ncpt;
    
    % Generate parses
    initMP = generate_random_parses(I,lib,K,verbose);
    nl = numel(initMP);
    init_scores = zeros(nl,1);
    for i=1:nl
        init_scores(i) = scoreMP_NoRel(initMP{i},lib);
    end      
    
    % Visualize initial parses
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

    % Run search for each candidate
    finalMP = cell(nl,1);
    for i=1:nl
        finalMP{i} = SearchForParse(initMP{i},lib,verbose);
    end
    final_scores = zeros(nl,1);
    for i=1:nl
       final_scores(i) = scoreMP(finalMP{i},lib); 
    end
    
    % Visualize the final parses
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
    
    % return structure
    G = struct;
    G.models = finalMP;
    G.scores = final_scores;
    G.init_scores = init_scores;
    G.init_models = initMP;
    G.img = I;
    G.PM = ps;
    
end