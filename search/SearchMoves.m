classdef SearchMoves < BetterHandle
    % SEACHMOVES... Operators for the search algorithm that finds
    %               the best-fitting motor program.
    %
    % It does not modify the coarse stroke decomposition, besides
    % by optimizing stroke direction and stroke split/merges
    %
    
    properties
        MAX_NS_ALL_PERM = 6; % consider all flips/permutations for characters
            % up to this size (-2 for fast mode)
    end
    
    properties
        M                
        searchPM
    end
    
    properties (Dependent)
        lib
        verbose
        fast_mode
    end
    
    methods
        
        function this = SearchMoves(Minit,lib,verbose,fast_mode)
            % constructor
            % 
            % Input:
            %  Minit : model to optimize
            %  lib: library
            %  verbose: (true/false) describe proress?
            %  fast_mode: (true/false) if true, do not run gradient-based optimization
            % 
            if ~exist('verbose','var')
               verbose = false; 
            end
            if ~exist('fast_mode','var')
                fast_mode = false; 
            end            
            assert(~Minit.has_relations);
            if fast_mode
                this.MAX_NS_ALL_PERM = this.MAX_NS_ALL_PERM - 2;
            end
            this.M = Minit;
            this.searchPM.lib = lib;
            this.searchPM.verbose = verbose;
            this.searchPM.fast_mode = fast_mode;
            this.searchPM.MAX_NS_ALL_PERM = this.MAX_NS_ALL_PERM;            
        end
        
        function out = get.lib(this)
            out = this.searchPM.lib;
        end
        
        function out = get.verbose(this)
            out = this.searchPM.verbose; 
        end
        
        function out = get.fast_mode(this)
           out = this.searchPM.fast_mode; 
        end
        
        function disp_score(this)
            % display the current model score
            if this.M.has_relations()
                score = scoreMP(this.M,this.lib);
            else
                score = scoreMP_NoRel(this.M,this.lib);
            end
            if this.verbose, fprintf(1,' score = %s\n',num2str(score,4)); end
        end
        
        function move_opt_subids(this,list_sid)
           % optimize all of the sub-ids in "list_sid" (default is all 1:M.ns)
           if exist('list_sid','var')
               optimize_subids(this.searchPM,this.M,list_sid); 
           else
               optimize_subids(this.searchPM,this.M); 
           end
        end
        
        function move_opt_grad(this,list_sid)
           % gradient optimization, for the strokes listed in "list_sid"
           if ~exist('list_sid','var')
               optimize_grad(this.searchPM,this.M); 
           else
               optimize_grad(this.searchPM,this.M,list_sid); 
           end
        end
        
        function move_opt_order(this)
           % optimize the stroke order
           optimize_order(this.searchPM,this.M);
        end
        
        function bool_flip = move_opt_direction(this,sid)
           % optimize the direction of a stroke
           assert(~this.M.has_relations());
           Q = this.M.copy();
           flip_direction(this.searchPM,Q,sid);
           scoreM = scoreMP_NoRel(this.M,this.lib);
           scoreQ = scoreMP_NoRel(Q,this.lib);
           bool_flip = false;
           if scoreQ > scoreM
              this.M = Q;
              bool_flip = true;
           end
        end
        
        function move_split_merge(this)            
           % run sub-search that tries split/merge moves 
           if ~this.fast_mode
                this.M = SearchSplitMerge(this.M,this.lib,this.verbose);
           end
        end
        
        function move_opt_dir_order_rel(this)
            % optimize the direction, order, and relations between strokes
          
            % flip each of the strokes
            if this.verbose, fprintf(1,'  try flipping each stroke direction.\n'); end
            Q_flip = cell(this.M.ns,1);
            for sid=1:this.M.ns
               Q_flip{sid} = this.M.copy();
               flip_direction(this.searchPM,Q_flip{sid},sid);
            end
            
            % try all combinations of flips, 
            % with the optimal stroke order for each
            if this.verbose, fprintf(1,'  find optimal directions/orders.\n'); end
            if this.M.ns <= this.MAX_NS_ALL_PERM
                bin = all_binary_strings(this.M.ns);
            else
                bin = rand(2^this.MAX_NS_ALL_PERM,this.M.ns)>.5;
                bin = unique(bin,'rows');
            end
            nb = size(bin,1);    
            scores = zeros(nb,1);
            store_Q = cell(nb,1);
            for i=1:nb
               flip = bin(i,:);
               Q = this.M.copy();
               for sid=1:this.M.ns
                  if flip(sid) % if we flipped this stroke
                    Q.S{sid} = Q_flip{sid}.S{sid}.copy();
                  end
               end               
               optimize_order(this.searchPM,Q);
               optimize_relations(this.searchPM,Q);
               scores(i) = scoreMP(Q,this.lib);
               store_Q{i} = Q;               
            end
            
            % select the best combination of direction flips/stroke order
            [~,windx] = randmax(scores);
            % windx = argmax(scores);
            this.M = store_Q{windx}.copy();            
        end
 
    end
    
    
    
end

function base = all_binary_strings(n)
% Generate all binary sequences of length n
% base is matrix, where rows are sequences
   base = [true; false];
   for i=2:n
      [nb,ncol] = size(base);
      base = [true(nb,1) base; false(nb,1) base];
   end
end

function ll = optimize_relations(searchPM,Q)
% plug in the optimal relations between the strokes
%
% ll : best log-likelihood during relation search, which includes
%      all CPDs that directly depend on the relations
    llvec = argmax_relations(searchPM.lib,Q);
    ll = sum(llvec);
end

function optimize_grad(searchPM,Q,list_sid)
% run gradient-based optimization, and return
% 
    if ~exist('list_sid','var')
       list_sid = 1:Q.ns; 
    end
    
    % run gradient-based optimization
    if ~searchPM.fast_mode
        scoreF = argmax_fit_type(Q,searchPM.lib,list_sid,searchPM.verbose);
    end
    
    % Check to make sure the object was properly updated
    if Q.has_relations
        scoreQ = scoreMP(Q,searchPM.lib,'strokes',list_sid,'type',true,'token',true,'image',true);
    else
        scoreQ = scoreMP_NoRel(Q,searchPM.lib,'strokes',list_sid,'type',true,'token',true,'image',true);
    end
    if ~searchPM.fast_mode
        assert( isinf(scoreF) || aeq(scoreF,scoreQ) );
    end
    
end

function optimize_order(searchPM,Q)
    % optimize the stroke order, where we implicitly maximize
    % over relations. But the function returns a Q where the relations are
    % not set

    if isfield(Q.S{1},'R')
       error('cannot optimize order after relations are set'); 
    end

    % get all the permutations to try
    if Q.ns <= searchPM.MAX_NS_ALL_PERM
        % try all possible combinations
        P = perms(1:Q.ns);
    else
        % try a subset of the permutations
        np = factorial(searchPM.MAX_NS_ALL_PERM); %720        
        P = zeros(np,Q.ns);
        for i=1:np
            P(i,:) = randperm(Q.ns);
        end
        P = unique(P,'rows');
    end
        
    % score all the permutations    
    n = size(P,1);
    scores = zeros(n,1);
    for i=1:n
       QQ = Q.copy(); 
       perm = P(i,:);       
       QQ.S = QQ.S(perm);       
       scores(i) = optimize_relations(searchPM,QQ);
    end

    % pick the best order, and return Q without instantiating relations
    [~,windx] = randmax(scores);
    perm = P(windx,:);
    Q.S = Q.S(perm);
end

function optimize_subids(searchPM,Q,list_sid)
    % for all strokes in list_sid (default is all of them),
    % apply one iteration of coordinate ascent on the
    % sub-stroke ids
    if ~exist('list_sid','var')
       list_sid = 1:Q.ns; 
    end
    for sid=list_sid % each stroke
       if searchPM.verbose, fprintf(1,'   choose subid for stroke %d ',sid); end
       optimize_this_subid(Q,sid,searchPM.lib,searchPM.verbose);
    end    
end

function flip_direction(searchPM,Q,sid)
% reverse the direction of a stroke "sid",
% an optimize various features about that stroke
    UtilMP.flip_stroke(Q.S{sid});
    optimize_subids(searchPM,Q,sid);
    optimize_grad(searchPM,Q,sid);
end