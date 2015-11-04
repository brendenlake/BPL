classdef SearchSubStk < BetterHandle
    %
    % SearchSubStk : Search moves to parse a trajectory into sub-strokes 
    %
    % Trajectory is divided in various places into sub-strokes.
    %   Each sub-stroke is fit with a spline, and then
    %   this spline is categorized and scored with the prior on strokes
    %
    % Constraint: the spline must be within a certain distance of 
    %   each point it is trying to fit. Thus, we cannot wander too
    %   far off the image skeleton
    %
    % Inherent tradeoff:
    %   Thus, two competing forces: the more partitions we use,
    %   the more likley we get a good fit with our library (and fit constraints).
    %   However, the complexity term from the prior keeps it from growing too complex.
    %
    properties
       hyp % current hypothesis (boolean, where true denotes splits)
    end
    
    properties (SetAccess = private)
        traj % original trajectory
        lib
        hist_hyp
        ps % parameters
        verbose
    end
    
    properties (Dependent = true)
        tlen % length of trajectory
        parse % [ns x 1 cell] trajectory broken into pieces
        parse_indx % [ns x 1 cell] indx (into traj) for elements in parse
        curr_score % score of current hypothesis
        feasible % is the current hypothesis feasible?
        curr_prob_bid % current sub-stroke that has larget error
    end
       
    methods
        
        % constructor
        %
        % traj: trajectory in image space
        function this = SearchSubStk(traj,lib,verbose)
            this.traj = traj;
            this.lib = lib;            
            this.hyp = false(1,this.tlen);
            this.hyp(1) = true;
            this.hyp(end) = true;
            this.ps = defaultps_bottomup;
            
            if ~exist('verbose','var')
               verbose = false; 
            end
            this.verbose = verbose;
        end
        
        % get length of trajectory
        function y = get.tlen(this)
            y = size(this.traj,1);
        end
        
        % get the current parse into sub-strokes
        function S = get.parse(this)
           S = this.make_parse(this.hyp); 
        end
        
        % get the corresponding indices for that parse
        function IDX = get.parse_indx(this)
           [~,IDX] = this.make_parse(this.hyp); 
        end
        
        % get the current score of the hypothesis
        function score = get.curr_score(this)
           score = this.score(this.hyp);            
        end
        
        function problem_bid = get.curr_prob_bid(this)
           [~,out] = this.helper_score(this.hyp);
           problem_bid = out.problem_bid;
        end
        
        % get a feasible solution 
        function y = get.feasible(this)
            y = ~isinf(this.curr_score);
        end
        
        % run the search procedure
        %
        % Output
        %  S: [nsub x 1 cell] parsed trajectory
        function S = run(this)
            
            if this.tlen == 1
               S = {this.traj};
               return
            end
                       
            this.find_feasible_hyp();
            if ~this.feasible % just return if failed to find feasible
                S = {this.traj};
                return
            end
            assert(this.feasible);
            this.hist_hyp = this.hyp;
            this.optimize_direction();
                        
            % loop until we can't improve the score
            move_count = 1;
            while true
                
                % try all moves
                prop_split = this.propose_splits(this.ps.ntry_split);
                prop_merge = this.propose_merges();
                prop_wiggle = this.propose_wiggles();
                prop_replace = this.propose_replace();
                prop_hyp = logical([prop_split; prop_replace; prop_merge; prop_wiggle]);
                
                % remove hypotheses we have tried
                prop_hyp = this.rmv_tried_hyp(prop_hyp);
                
                ll = this.score(prop_hyp);
                
                % visualize if verbose
                if this.verbose
                   fprintf(1,'move %d\n',move_count);                   
                   fprintf(1,'splits.\n');
                   this.viz_candidates(prop_split);
                   fprintf(1,'replace.\n');
                   this.viz_candidates(prop_replace);
                   fprintf(1,'merges.\n');
                   this.viz_candidates(prop_merge);
                   fprintf(1,'wiggles.\n');
                   this.viz_candidates(prop_wiggle);
                end                
                
                % select the best move
                [best_score,windx] = randmax(ll);
                if best_score > this.curr_score + 1 % accept
                    this.hyp = prop_hyp(windx,:);
                    this.hist_hyp = [this.hist_hyp; this.hyp];
                    this.optimize_direction();
                    move_count = move_count + 1;
                else % no improvement
                    if this.verbose
                       fprintf(1,'search complete\n');
                       this.viz_search(); 
                    end
                    break; 
                end
                
            end
            
            % return the best parse
            S = this.parse;
           
        end
        
        % flip the stroke direction
        function flip_dir(this)
            this.hyp = this.hyp(end:-1:1);
            this.hist_hyp = this.hist_hyp(:,end:-1:1);
            this.traj = this.traj(end:-1:1,:);
        end
        
        % optimize the direction of the stroke
        function optimize_direction(this)
            curr_score = this.curr_score;
            this.flip_dir;
            flip_score = this.curr_score;
            if aeq(flip_score,curr_score) % if scores are equal
                if rand > 0.5 % flip randomly
                    this.flip_dir;
                end                
            elseif flip_score > curr_score
                % do nothing if we like the current score
            else
                this.flip_dir;
            end
        end
        
        %
        % get the parse of the grand trajectory into
        % sub-stroke trjacetories
        %
        % Output
        %  S: [ns x 1 cell] parse into trajectories
        %  IDX: [ns x 1 cell] same as S, except it contains
        %   indices into the original parse
        function [S,IDX] = make_parse(this,myhyp)
            if numel(myhyp)==1
               S ={this.traj};
               IDX = {1};
               return
            end            
            assert(islogical(myhyp));
            b = find(myhyp);
            nsub = numel(b)-1;
            S = cell(nsub,1);
            IDX = cell(nsub,1);
            for bindx=1:nsub
                istart = b(bindx);
                iend = b(bindx+1);
                S{bindx} = this.traj(istart:iend,:);
                IDX{bindx} = vec(istart:iend);
            end
        end
        
        %
        % propose merging two sub-strokes in 
        % every possible place
        %
        % Output
        %  prop_hyp: [m x tlen boolean] rows are hypotheses
        function prop_hyp = propose_merges(this)
            b = find(this.hyp);
            b(1) = [];
            b(end) = [];
            nb = numel(b);            
            prop_hyp = [];
            for bid=1:nb
               tormv = b(bid); 
               prop_hyp(bid,:) = this.hyp;
               prop_hyp(bid,tormv) = false;
            end
            prop_hyp = logical(prop_hyp);
        end
        
        %
        % sample a set of possible splits for a trajectory
        %
        % Output
        %  prop_hyp: [nsamp x tlen boolean] rows are hypotheses
        %     may have fewer rows if duplicates are chosen
        %
        function prop_hyp = propose_splits(this,nsamp)            
           prop_hyp = false(nsamp,this.tlen);
           for i=1:nsamp
               prop_hyp(i,:) = this.helper_propose_split();
           end 
           prop_hyp = unique(prop_hyp,'rows');
        end
        
        %
        % sample new positions for the trajectory breaks
        %
        % Output
        %  prop_hyp: [nsamp x tlen boolean] rows are hypotheses
        %     may have fewer rows if duplicates are chosen
        %  nsamp: number of samples we want to take
        function prop_hyp = propose_wiggles(this)
            b = find(this.hyp);
            b(1) = [];
            b(end) = [];
            nb = numel(b);            
            prop_hyp = [];
            for bid=1:nb
               prop_hyp(bid,:) = this.helper_propose_wiggle(bid); 
            end 
            prop_hyp = logical(prop_hyp);
        end
        
        %
        % Propose replacing each sub-stroke break
        % by removing it and then re-sampling it
        %
        function prop_hyp = propose_replace(this)
            b = find(this.hyp);
            b(1) = [];
            b(end) = [];
            nb = numel(b);            
            prop_hyp = [];
            for bid=1:nb
               prop_hyp(bid,:) = this.helper_propose_replace(bid); 
            end 
            prop_hyp = logical(prop_hyp);
        end        
        
        %
        % score the hypothesis
        %
        % Input
        %  list_hyp: [n x tlen] rows are different hypotheses to score
        %
        % Output
        %  ll [n x 1]: approximation of log-prior for this stroke
        %  cell_out [n x 1 cell]: misc. output information
        %
        function [ll,cell_out] = score(this,list_hyp)
           
            % check inputs
            if isvector(list_hyp)
                assert(numel(list_hyp)==this.tlen);
                list_hyp = list_hyp(:)';
            else
                % rows must be hypotheses
                assert(size(list_hyp,2)==this.tlen); 
            end
            
            % score each one
            n = size(list_hyp,1);
            ll = zeros(n,1);
            cell_out = cell(n,1);
            for i=1:n
                [ll(i),cell_out{i}] = this.helper_score(list_hyp(i,:));
            end
            
        end
        
        %
        % visualize the search process
        %
        function viz_search(this)
            nstep = size(this.hist_hyp,1);
            figure(110);
            clf
            nrow = ceil(sqrt(nstep));            
            ll = zeros(nstep,1);
            cell_out = cell(nstep,1);
          
            % plot the candidates
            for i=1:nstep
               step_hyp = this.hist_hyp(i,:);
               [ll(i),cell_out{i}] = this.helper_score(step_hyp);
               subplot(nrow,nrow,i);
               viz_parse(cell_out{i}.hyp_parse,cell_out{i}.recon_parse);
               nsub = numel(cell_out{i}.hyp_parse);
               lab = ['score = ',num2str(ll(i),3),', error = ',num2str(cell_out{i}.max_error,3)];
               xlabel(lab);
               title(['step ',num2str(i),', nsub=',num2str(nsub)]);
            end
            
            figure(111)
            clf;
            plot(1:nstep,ll);
            ylabel('score');
            xlabel('steps');            
        end
        
    end
    
    methods (Access = private)
       
       % remove hypotheses we have already troed 
       function prop_hyp = rmv_tried_hyp(this,prop_hyp)  
           assert(size(prop_hyp,2)==this.tlen);
           [C,tormv] = intersect(prop_hyp,this.hist_hyp,'rows');
           if ~isempty(tormv)
              z=1; 
           end
           prop_hyp(tormv,:) = [];         
       end        
        
       % visualize the candidate hypotheses 
       function viz_candidates(this,prop_hyp)
            nstep = size(prop_hyp,1);
            nshow = nstep + 1;
            figure(112);
            clf
            nrow = ceil(sqrt(nshow));
            
            % plot the current hypothesis
            subplot(nrow,nrow,1);
            [ll_curr,out_curr] = this.helper_score(this.hyp);
            viz_parse(out_curr.hyp_parse,out_curr.recon_parse);
            lab = ['score = ',num2str(ll_curr,3),', error = ',num2str(out_curr.max_error,3)];
            xlabel(lab);
            nsub = numel(out_curr.hyp_parse);
            title(['Current, nsub=',num2str(nsub)]);
            
            % plot the moves
            ll = zeros(nstep,1);
            cell_out = cell(nstep,1);
            for i=1:nstep
               step_hyp = prop_hyp(i,:);
               [ll(i),cell_out{i}] = this.helper_score(step_hyp);
               subplot(nrow,nrow,i+1);
               viz_parse(cell_out{i}.hyp_parse,cell_out{i}.recon_parse);
               nsub = numel(cell_out{i}.hyp_parse);
               lab = ['score = ',num2str(ll(i),3),', error = ',num2str(cell_out{i}.max_error,3)];
               xlabel(lab);
               title(['option ',num2str(i),', nsub=',num2str(nsub)]);
            end
            input('press enter to continue\n','s');
       end
        
       % Find a feasible hypothesis to begin with 
       function find_feasible_hyp(this)
           count = 0;
           while ~this.feasible             
               this.hyp = this.helper_propose_split(this.curr_prob_bid);
               count = count + 1;
               if count > 100
                   return;
                   %error('unable to find a feasible starting location');
               end
           end           
       end    
        
       % reconstruct each segment of the parse 
       function S = reconstruct_parse(this,stats)
            nsub = numel(stats.ids);
            S = cell(nsub,1);
            hyp_parse = this.make_parse(stats.hyp);
            for bid=1:nsub
                
                % correct scale
                X = stats.shapes(:,:,bid) .* stats.invscales(bid);
                
                % correct position
                inv_com = -stats.com(bid,:);
                X = offset_stk(X,inv_com);
                
                % evaluate spline
                neval = size(hyp_parse{bid},1);
                S{bid} = get_stk_from_bspline(X,neval);
            end
       end
       
       
       % propose wiggling one of the current sub-stroke break
       % points
       %
       % bindx: (optional) select which sub-stroke break to wiggle,
       %   or choose randomly if absent
       function prop_hyp = helper_propose_wiggle(this,bindx)
          
            % sample the divider we want to wiggle
            b = find(this.hyp);
            b(1) = [];
            b(end) = [];
            nb = numel(b);
            assert(nb >= 1);
            if ~exist('bindx','var')
                bindx = randint(1,1,[1 nb]);
            else
                assert(bindx <= nb);
            end
            
            % sample the new location of the break point
            curr_loc = b(bindx);
            new_loc = curr_loc + this.sample_wiggle;
            prop_hyp = this.hyp;
            if new_loc >= 1 && new_loc <= this.tlen
                prop_hyp(curr_loc) = false;
                prop_hyp(new_loc) = true;
            end
            
       end
       
       % sample an offset for some hypothetical
       % sub-stroke break
       function shift = sample_wiggle(this)
           sd = round(this.ps.sigma_wiggle);
           x = [-3*sd:-1, 1:3*sd];
           logy = mvnormpdfln(x,0,sd);
           logy = logy - logsumexp(logy(:));
           py = exp(logy);
           py = py ./ sum(py);
           indx = find(mnrnd(1,py));
           assert(numel(indx)==1);
           shift = x(indx);
       end
       
       % propose replacing a sub-stroke break with another one
       % 
       % Input
       %  bindx: (optional) the precise break we want to re-split
       function prop_hyp =  helper_propose_replace(this,bindx)
            
           % sample the divider we want to replace
            b = find(this.hyp);
            b(1) = [];
            b(end) = [];
            nb = numel(b);
            assert(nb >= 1);
            if ~exist('bindx','var')
                bindx = randint(1,1,[1 nb]);
            else
                assert(bindx <= nb);
            end
           
            prop_hyp = this.hyp;
            curr_loc = b(bindx);
            prop_hyp(curr_loc) = false;
           
            
            [S,IDX] = this.make_parse(prop_hyp);
            nsub = numel(S);
            
            % check to make sure sub-stroke is big enough
            buffer = round(this.lib.ncpt./2);
            subtraj = S{bindx};
            len = size(subtraj,1);
            if len < 2*buffer
               return; % we can't do anything, so just return 
            end
            
            % choose split index
            isplit = sample_split(subtraj,buffer);
            itotal = IDX{bindx}(isplit);
            prop_hyp(itotal) = true;
           
       end
       
       
       % propose splitting a sub-stroke
       %
       % Input
       %   bid: (optional) the precise sub-stroke we should try splitting
       function prop_hyp = helper_propose_split(this,bid)
            prop_hyp = this.hyp;
            [S,IDX] = this.make_parse(this.hyp);
            nsub = numel(S);
            if exist('bid','var')
                assert(bid <= nsub);
            else
                bid = randint(1,1,[1 nsub]);
            end
            
            % check to make sure sub-stroke is big enough
            buffer = round(this.lib.ncpt./2);
            subtraj = S{bid};
            len = size(subtraj,1);
            if len < 2*buffer
               return; % we can't do anything, so just return 
            end

            % choose split index
            isplit = sample_split(subtraj,buffer);
            itotal = IDX{bid}(isplit);
            prop_hyp(itotal) = true;
        end
       
        % score the hypothesis
        %
        % Output
        %  score: approximation of log-prior for this stroke,
        %     where it is -inf if we violate closeness constraints
        %  out.max_error: maximum point-wise error across reconstructions
        %
        function [score,out] = helper_score(this,hyp)            
            assert(isvector(hyp));
            assert(islogical(hyp));
            assert(hyp(1)); % must be bounded at ends
            assert(hyp(end));
            hyp = hyp(:);
            S = this.make_parse(hyp);
            
            % cluster the sub-strokes as primitives
            nsub = length(S);
            [cell_shapes,cell_com,cell_scales] = ...
                normalize_dataset(S,this.lib.newscale,false);
            clusterid = cell(nsub,1);
            for b=1:nsub
                clusterid{b} = classify_traj_as_subid(cell_shapes{b},cell_scales{b},this.lib);
            end
            
            % extract their shape, scales, and score
            stats = struct;
            stats.hyp = hyp;
            stats.ids = zeros(nsub,1);
            stats.invscales = zeros(nsub,1);
            for b=1:nsub
              stats.ids(b) = clusterid{b}.indx;
              stats.invscales(b) = 1 ./ cell_scales{b}(1);
              stats.shapes(:,:,b) = clusterid{b}.bspline;
              stats.com(b,:) = cell_com{b};
            end
            ll = this.score_stroke(stats);
           
            % Compute erorr on reconstruction of hypothesis
            hyp_parse = this.make_parse(hyp); % parse real trajectory
            recon_parse = this.reconstruct_parse(stats);            
            [max_error,problem_bid] = compute_parse_dist(hyp_parse,recon_parse);
           
            % output structure
            out = struct;
            out.ll = ll;
            out.max_error = max_error;
            out.problem_bid = problem_bid;
            out.hyp_parse = hyp_parse;
            out.recon_parse = recon_parse;
        
            % see if we follow the constraints
            score = ll;
            if max_error > this.ps.traj_abs_error_lim
                score = -inf;
            end
            
        end
        
        % Score a stroke with the prior objective function
        function ll = score_stroke(this,stats)
           ll = 0;
           ll = ll + CPD.score_sequence(this.lib,-1,stats.ids);
           ll = ll + sum(CPD.score_invscale_type(this.lib,stats.invscales,stats.ids));
           ll = ll + sum(CPD.score_invscale_token(this.lib,stats.invscales,stats.invscales));
           ll = ll + sum(CPD.score_shape_marginalize(this.lib,stats.shapes,stats.ids));
        end
       
    end
    
end

% compare parse of original trajcetory
% to the spline reconstruction
function viz_parse(parse,parse_recon)
    hold on
    
    nsub = numel(parse);
    assert(nsub == numel(parse_recon));

    for bid=1:nsub
        seg = parse{bid};
        plot(seg(:,1),seg(:,2),'k','LineWidth',2);
        plot(seg(end,1),seg(end,2),'.k','MarkerSize',24);
    end
    
    for bid=1:nsub
        seg = parse_recon{bid};
        plot(seg(:,1),seg(:,2),'r','LineWidth',1);
        plot(seg(end,1),seg(end,2),'.r','MarkerSize',16);
    end
   
    
end

% find, over all points, the maximum distance for which a spline approximation
% is off.
%
% Input
%  cell_t1: [nsub x 1 cell] parse of trajectories
%  cell_t2: 
%
% Output
%  grand_max: [scalar] answer
%  problem_b: [scalar] which cell does the problem point come from
%
function [grand_max,problem_bid] = compute_parse_dist(cell_t1,cell_t2)
    nsub = numel(cell_t1);
    assert(numel(cell_t2)==nsub);
    max_dist = zeros(nsub,1);
    for b=1:nsub
       max_dist(b) = compute_dist(cell_t1{b},cell_t2{b});        
    end     
    [grand_max,problem_bid] = max(max_dist);
end

% average distance between trajectories
function [max_dist,ave_dist] = compute_dist(t1,t2)
    n1 = size(t1,1);
    n2 = size(t2,1);
    assert(n1==n2)
    dist = zeros(n1,1);
    for i=1:n1
       v1 = t1(i,:);
       v2 = t2(i,:);
       dist(i) = norm(v1(:)-v2(:));
    end
    ave_dist = mean(dist);
    max_dist = max(dist);
end

% Sample a split point, where the probability
% is proportional to the square of the angle
function isplit = sample_split(motor,buffer)
    n = size(motor,1);
    theta = zeros(n,1);
    
    % compute the angle between successive points
    for i=3:n-2
        v1 = motor(i,:) - motor(i-2,:);
        v2 = motor(i+2,:) - motor(i,:);
        cos_theta = (v1(:)' * v2(:)) ./ (max(norm(v1),1) * max(norm(v2),1)) ;
        theta(i) = acosd(cos_theta);
    end    
    pvec = theta.^2;
    pvec = pvec + 1e-4;
    pvec = pvec ./ sum(pvec);
    
    % zero out the buffer
    pvec(1:buffer-1) = 0;
    pvec(end-buffer+1:end) = 0;
    pvec = pvec ./ sum(pvec);
    assert(~any(isnan(pvec)));    
    
    % select point to split on
    isplit = find(mnrnd(1,pvec));    
end