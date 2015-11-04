classdef Walker < BetterHandle
    % Walker: Class that controls a walk on the directed graph
    %   that is defined by the skeleton. The walk is complete
    %   when all of the edges have been covered.
    %
    
    properties
       list_WS = []; % list of WalkerStroke objects   
       T % undirected graph structure
    end
    
    properties (Dependent = true)
       ns % number of strokes
       S % full trajectories for each stroke
       edges_visited % [ne x 1 bool] is each edge visited?
       nodes_visited % [n x 1 bool] is each node visited?
       complete % [scalar] is the entire graph drawn?
       curr_ni  % current nodes index
       curr_pt  % current point location
    end
    
    methods
        
        % 
        function this = Walker(T)
           this.T = T;
        end
       
        % get the number of strokes
        function y = get.ns(this)
           y = numel(this.list_WS); 
        end
        
        % if all of the edges are visited, 
        % we are done drawing
        function y = get.complete(this)
           y = all(this.edges_visited); 
        end
        
        % get boolean list of all edges visited
        function edges = get.edges_visited(this)
           ne = numel(this.T.S);
           edges = false(ne,1);
           for sid=1:this.ns
               WS = this.list_WS{sid};
               edges(WS.visited_ei) = true;
           end
        end
        
        % get boolean list of all nodes visited
        function nodes = get.nodes_visited(this)
           n = this.T.n;
           nodes = false(n,1);
           for sid=1:this.ns
              WS = this.list_WS{sid};
              nodes(WS.visited_ni) = true; 
           end
        end
        
        % get current node index
        function ni = get.curr_ni(this)
           WS = this.list_WS{end}; 
           ni = WS.curr_ni; 
        end
       
        % get the current node point
        function pt = get.curr_pt(this)
           pt = this.T.G(this.curr_ni,:);
        end
        
        % Add all singleton nodes in the skeleton graph,
        % that have no edges.
        function add_singletons(this)
            for i=1:this.T.n % for each node
                if sum(this.T.E(i,:))==0 % if singleton
                   start_pt = this.T.G(i,:);
                   this.list_WS{end+1} = WalkerStroke(this.T);
                   this.list_WS{end}.start_pt = start_pt;
                end
            end
        end
         
        % get the full trajectory of the parse
        %
        % Output
        %  S: [ns x 1 cell] of parse
        function S = get.S(this)
           S = cell(this.ns,1);
           for sid=1:this.ns
               S{sid} = this.list_WS{sid}.stk;
           end
        end
        
        % Plot the current walk
        function plot(this)
           myplot(this.S); 
        end
        
        % Accessor for plot_traj
        function h = plot_path(this,path,color)
           h = plot_traj(path,color);
        end

        % Get all the possible local moves we can make.
        function [cell_traj,vei] = get_moves(this)
            [cell_traj,vei] = this.list_WS{end}.get_moves();
        end
        
        % Select one of the moves
        function select_moves(this,sel)
           this.list_WS{end}.select_move(sel);
        end
        
        % Get local moves we can make that cover a new edge
        %
        % Input
        %  pt: [1x2] optional, compute moves from this point
        function [cell_traj,vei] = get_new_moves(this)
            [cell_traj,vei] = this.list_WS{end}.get_moves();
            n = numel(vei);
            
            % select items that were not visited
            sel = false(n,1);
            for i=1:n
               if ~this.edges_visited(vei(i))
                  sel(i) = true; 
               end
            end
            
            vei = vei(sel);
            cell_traj = cell_traj(sel);
        end
        
        % Select local moves we can make that cover a new edge
        function select_new_moves(this,sel_new)
            [~,vei] = this.list_WS{end}.get_moves();
            
            % select items that were not visited
            n = numel(vei);
            mysel = false(n,1);
            for i=1:n
               if ~this.edges_visited(vei(i))
                  mysel(i) = true; 
               end
            end
            
            findx = find(mysel);
            sel = findx(sel_new);
            this.list_WS{end}.select_move(sel);
        end
        
    end
   
end

% Plot trajectory
function myplot(S)
    hold on
    ns = numel(S);
    for i=1:ns
        stk =S{i};
        col = get_color(i);
        plot_traj(stk,col);
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
    h = plot(stk(:,1),stk(:,2),col,'LineWidth',4);
end

% Color map for the stroke of index k
function out = get_color(k)
    scol = {'r','g','b','m','c'};
    ncol = length(scol);
    if k <=ncol
       out = scol{k}; 
    else
       out = scol{end}; 
    end
end