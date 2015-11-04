classdef RandomWalker < Walker
    %
    % RandomWalker : produce a random walk on the graph skeleton.
    %  It is presumed the direction and order of the strokes doesn't
    %  matter, and this will be optimized later
    %
    %  Methods
    %   .sample(nsamp) : number of samples you want to make
    %   .det_walk      : determinisitc walk following minimal angle
    %   .make(true)    : walks you through a walk, visually. for debugging.
    %
    % Steps for generating a random walk
    %   Step 1) Add singleton nodes in the graph
    % Repeat until complete
    %   Step 2) For each point on a new edge, place pen
    %      down on that point with weight inversely related to the 
    %      degree of new edges coming from it.
    %   Step 3) Choose a new node based on the curvature of the local spline fit.
    %
    properties
        ps
        verbose
        exp_wt_start
        lambda_softmax
    end
    
    methods
        
        % Constructor
        %
        % Input
        %  T: undirected skeleton graph
        function this = RandomWalker(T)
            % superclass constructor
            this = this@Walker(T);
            this.ps = defaultps_bottomup;
        end        
        
        % Produce "nsamp" samples
        % from the random walk model
        function list_samples = sample(this,nsamp)
            list_samples = cell(nsamp,1);
            for i=1:nsamp
               myns = inf;
               
               % sample hyper-parameters
               this.exp_wt_start = rand_discrete(this.ps.int_exp_wt);
               this.lambda_softmax = rand_discrete(this.ps.int_lambda_soft);
               
               while myns > this.ps.max_len
                   samp = this.make();
                   myns = numel(samp);
               end
               list_samples{i} = samp;
            end            
        end
        
        % Produce a deterministic random walk
        function S_walk = det_walk(this)
            
            % sample hyper-parameters
            this.exp_wt_start = 1000;
            this.lambda_softmax = 1000;
            
            S_walk = this.make;
        end
        
        % Make a random walk through the graph
        function S_walk = make(this,verbose)
           if ~exist('verbose','var') 
              verbose = false;              
           end
           this.verbose = verbose;
           if verbose, h = figure; end
           
           % plot skeleton
           if verbose, this.T.plot_skel(); end           
           this.clear(); % clear all work
           this.add_singletons();
           if ~this.complete    
               this.pen_up_down();
               if verbose
                   this.plot();
                   input('press <enter> to continue\n','s');
               end
           end
           
           % continue until complete
           while ~this.complete
              this.pen_angle_step();
              if verbose
                  figure(h);
                  this.plot();
                  input('press <enter> to continue\n','s');
              end
           end
           
           S_walk = this.S;
        end
        
        % Clear the object
        function clear(this)
           this.list_WS = []; 
        end               
                        
        % Place you pen down at an unvisited edge,
        % inversely proportional to the number of unvisited
        % paths going from it.
        function pen_up_down(this)
            [new_pts,degree] = this.pts_on_new_edges;
            logwts = this.exp_wt_start .* log(1./degree);
            logwts = logwts - logsumexp(logwts(:));
            wts = exp(logwts);
            [~,rindx] = rand_discrete(wts,wts);
            start_pt = new_pts(rindx,:);            
            this.list_WS{end+1} = WalkerStroke(this.T);
            this.list_WS{end}.start_pt = start_pt;
            if ~this.complete
               this.pen_simple_step(); 
            end
        end
            
        % Angle move: selet a step based on the 
        %   angle from the current trajectory
        function pen_angle_step(this)            
            [cell_traj,vei] = this.get_moves();
            if isempty(vei)
               this.pen_up_down();
               return
            end

            % select edges that were not visited yet
            n = numel(vei);
            newedge = false(n,1);
            for i=1:n
               if ~this.edges_visited(vei(i))
                  newedge(i) = true; 
               end
            end
            
            % Default angle for used edges
            angles = this.ps.faux_angle_repeat * ones(n,1);
            [angles_for_new,first_half,second_half,cell_smooth] = angles_for_moves(this,cell_traj(newedge));
            angles(newedge) = angles_for_new;
            angles = [angles; this.ps.faux_angle_lift];
            
            % Select move stochastically
            rindx = action_via_angle(this,angles);
            if rindx == numel(angles) % if we chose to lift
                this.pen_up_down();
            else % regular step
                this.select_moves(rindx); 
            end  
            
            % Visualize the move if desired
            if this.verbose
                angles
                this.viz_angle_step(angles_for_new,first_half,second_half,cell_smooth);
                pause(.1);
            end
            
        end
        
        % Visualize an angle-based move if we want to
        function viz_angle_step(this,angles_for_new,first_half,second_half,cell_smooth) 
            n = numel(angles_for_new);
            nrow = ceil(sqrt(n)); 
            figure(10);
            clf
            for i=1:n
                subplot(nrow,nrow,i);
                hold on
                this.T.plot_skel();
                first_T = first_half{i};
                second_T = second_half{i};
                smooth_T = cell_smooth{i};
                plot_traj(smooth_T,'y');
                plot_traj(first_T,'c');
                plot_traj(second_T,'m');
                
                for j=1:this.T.n
                    plot(this.T.G(j,2),this.T.G(j,1),'r.','MarkerSize',7); 
            	end
                
                angle = angles_for_new(i);
                xlabel(num2str(angle,3));
            end
        end
        
        % Simple move: select a step uniformly at random
        %   from the step of new edges. Do not consider lifting
        %   the pen until you run out of new edges.
        function pen_simple_step(this)
            [~,vei] = this.get_new_moves();
            if isempty(vei)
               this.pen_up_down();
               return
            end
            sel = randint(1,1,[1 numel(vei)]);
            this.select_new_moves(sel); 
        end
        
        
    end
    
    methods (Access = private)        
        
        % Action model given a vetor of "angles" (real or faux)
        %
        % Move probability proportional to exp(-lambda*angle/180)
        %
        function rindx = action_via_angle(this,angles)
           assert(isvector(angles));
           theta = angles(:) ./ 180;
           netinput = -this.lambda_softmax*theta;
           logpvec = netinput - logsumexp(netinput);
           pvec = exp(logpvec);        
           n = numel(angles);           
           [~,rindx] = rand_discrete(1:n,pvec);
        end

        % For each move, compute its angle.
        %
        % Input
        %  cell_traj: [k x 1 cell]
        %  
        % Output
        %  angles: [nt x 1 scalar] angles in degrees
        %  first_half: [nt x 1 cell]
        %  second_half: [nt x 1 cell]
        function [angles,first_half,second_half,cell_smooth] = angles_for_moves(this,cell_traj)
            junct_pt = this.curr_pt;

            % For each possible move, list the entire stroke
            % that we would create if we accepted it.
            last_stk = this.S{end};            
            nt = numel(cell_traj); 
            cell_prop = cell(nt,1);
            for i=1:nt               
                cell_prop{i} = [last_stk; cell_traj{i}(2:end,:)];
            end

            % Smooth each candidate stroke
            cell_smooth = cell(nt,1);
            for i=1:nt
               cell_smooth{i} = fit_smooth_stk(cell_prop{i},this.T.I,false,this.ps);
            end

            % At the junction, isolate the relevant segments of the
            % smoothed stroke
            first_half = cell(nt,1);
            second_half = cell(nt,1);
            angles = zeros(nt,1);
            for i=1:nt
               [first_half{i},second_half{i}] = split_by_junction(junct_pt,cell_smooth{i},this.ps.rad_junction);
               angles(i) = compute_angle(second_half{i},first_half{i},this.ps);
            end
            if any(isnan(angles) | imag(angles)~=0)
                error('error in angle calculation');
            end
        end  
                  
        % For all the new edges in the graph, make a list
        % of their start/end points where we may want to drop our pen.
        % Also, return their degree (counting only unvisited edges)
        %
        % Output
        %  list_pts: [n x 2] points
        %  degree:   [n x 1] degree (of new edges) from each of
        %            these points
        function [list_pts,degree] = pts_on_new_edges(this)
            new_edges = ~this.edges_visited; 
            fei = find(new_edges);
            nei = numel(fei);
            list_pts = [];
            for i=1:nei
               traj = this.T.S{fei(i)};
               list_pts = [list_pts; traj(1,:)];
               list_pts = [list_pts; traj(end,:)];
            end
            
            % compute the "new" degree for each option
            npt = size(list_pts,1);
            degree = zeros(npt,1);
            for i=1:npt
                vei = this.T.get_branches(list_pts(i,:));
                tormv = this.edges_visited(vei);                
                vei(tormv) = [];
                degree(i) = numel(vei);
            end   
        end 
        
    end
    
end

% plot a stroke trajectory in image space,
% where the x and y dimensions are reversed
function h = plot_traj(stk,col)
    ystk = stk(:,2);
    stk(:,2) = stk(:,1);
    stk(:,1) = ystk;
    if size(stk,1)==1
        h = plot(stk(1),stk(2),[col,'.'],'MarkerSize',24);
        return
    end   
    h = plot(stk(:,1),stk(:,2),col,'LineWidth',2);
end

% Get portion of trajectory within the specific radius,
% and divide it in two based on the closest point to the junction
%
% Input
%  junct_pt: [1 x 2] junction point in skeleton
%  traj:     [k x 2] smoothed trajectory
%  radius:   [scalar] radius around junction for computing angle
%
% Output
%  first_half_traj: [m x 2]
%  second_half_traj: [n x 2]
%
function [first_half_traj,second_half_traj] = split_by_junction(junct_pt,traj,radius)
    assert(numel(junct_pt)==2);
    assert(size(traj,2)==2);    
    d = pdist2(junct_pt,traj);
    
    % get the relevant trajectory portion
    sel = d < radius;
    
    % If we are re-tracing a previously visited junction, we
    % must be careful not to include any of the trajectory points
    % on this junction.
    last = find(sel,1,'last');
    for i=last-1:-1:1
       if ~sel(i+1)
           sel(i) = false;
       end
    end
    new_traj = traj(sel,:);
    
    % divide it into two halves
    bindx = argmin(d(sel));
    first_half_traj = new_traj(1:bindx,:);
    second_half_traj = new_traj(bindx:end,:);
    
end

% Compute the angle between two vectors.
%
% Output
%   theta: [scalar] between 0 and 180 degrees
function thetaD = compute_angle(seg_ext,seg_prev,ps)

    % make sure the vectors are long enough so we can
    % compute an angle
    n_ext = size(seg_ext,1);
    n_prev = size(seg_prev,1);
    if n_ext < 2 || n_prev < 2
       thetaD = ps.faux_angle_too_short;
       return
    end
    
    % compute angle in degrees
    v_ext = seg_ext(end,:) - seg_ext(1,:);
    v_prev = seg_prev(end,:) - seg_prev(1,:);
    v_ext = v_ext(:);
    v_prev = v_prev(:);
    denom = (norm(v_ext)*norm(v_prev));
    if aeq(denom,0)
       denom = 1; 
    end
    val = (v_ext'*v_prev) ./ denom;
    assert(isscalar(val));
    if val > 1, val = 1; end
    if val < -1, val = -1; end
    thetaD = acosd(val);
end