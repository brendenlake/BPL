classdef WalkerStroke < BetterHandle
   % WalkerStroke: A single stroke in a parse of a graph skeleton
    
   properties 
      start_pt  % [1 x 2] where the stroke starts (can be redundant)
                %         only necessary when ei is empty
      ei        % [k x 1] index vector of sequence of edges
      ei_flip   % [k x 1 bool] do we flip each edge, or leave as in T.S?
      T   % Undirected graph structure
   end
  
   properties (Dependent = true)
      stk % [n x 2] stroke corresponding to other variables
      k % number of edges
      curr_pt % current point on skeleton [1 x 2]
      curr_ni % index of that point [scalar]
      visited_ni % visited node indices [n x 1 boolean]
      visited_ei % visited edge indices [KK x 1 boolean] 
   end
   
   methods
       
        function this = WalkerStroke(T)
            this.T = T;
        end
        
        % number of steps
        function k = get.k(this)
           k = numel(this.ei); 
        end
        
        % get the point associated with the last node
        %
        % y: [1 x 2] point on plane
        function y = get.curr_pt(this)
           if isempty(this.ei)
              y = this.start_pt;
              return
           end            
           traj = this.T.S{this.ei(end)};
           if this.ei_flip(end)
              traj = flip(traj); 
           end
           y = traj(end,:);
        end
        
        % get the point associated with the last node
        %
        % ni: [scalar] index into T.G
        function ni = get.curr_ni(this)
            pt = this.curr_pt;
            ni = this.T.map_pts_to_ni(pt);
        end
       
        % which nodes did we visit?
        %
        % visited: [n x 1 boolean]
        function visited = get.visited_ni(this)
           visited = false(this.T.n,1);
           link_ei_to_ni = this.T.link_ei_to_ni; 
           array_ni = link_ei_to_ni(this.ei,:);
           vni = array_ni(:);
           visited(vni) = true;
        end
        
        % which edge did we visit?
        %
        % visited: [ne x 1 boolean]
        function visited = get.visited_ei(this)
           nei = numel(this.T.S);
           visited = false(nei,1); 
           visited(this.ei) = true; 
        end
        
        % Get the stroke trajectory in its full glory
        %
        % Output
        %  vstk: [n x 2] trajectory
        function vstk = get.stk(this)
           vstk = this.get_path();
        end
        
        % Given the current point, where can we go?
        %
        % Ouput
        %  cell_traj: [m x 1 cell] trajectory fragments we can choose from
        %  vei: [m x 1] list of edges we can go to
        function [cell_traj,vei] = get_moves(this)
           [vei,vei_flip] = this.T.get_branches(this.curr_pt);
           n = numel(vei);
           cell_traj = cell(n,1);
           for i=1:n
               vindx = vei(i);
               traj = this.T.S{vindx};
               if vei_flip(i)
                  traj = flip(traj); 
               end
               cell_traj{i} = traj;
           end           
        end
        
        % Select a move (sel) based on the list
        % returned from this.get_move
        %
        % Input
        %  sel: index into the list of moves from this.get_moves
        function select_move(this,sel)
            assert(isscalar(sel));
            [vei,vei_flip] = this.T.get_branches(this.curr_pt);
            nmove = numel(vei);
            assert(sel<= nmove && sel >= 1);
            this.ei = [this.ei; vei(sel)];
            this.ei_flip = [this.ei_flip; vei_flip(sel)];
        end
   
   end
   
   methods (Access = private)
              
        % Get the [x,y] path through the skeleton
        %
        % Output
        %  vstk: [nn x 2] path through the skel
        function vstk = get_path(this)
            
           if isempty(this.ei)
               vstk = this.start_pt;
               return
           end
           
           % first step
           traj = this.T.S{this.ei(1)};
           if this.ei_flip(1) % flip if needed
              traj = flip(traj);
           end
           vstk = traj;
            
           % for the additional steps
           for step=2:this.k
              traj = this.T.S{this.ei(step)};
              if this.ei_flip(step) % flip if needed
                  traj = flip(traj);
              end
              
              % make sure this is a valid step
              assert(isequal(vstk(end,:),traj(1,:)));
              
              % add the next step
              vstk = [vstk; traj(2:end,:)];
           end
           
        end
        
   end
    
end

% Flip direction of trajectory
%
% Input
%  traj: [n x 2] trajectory
function traj = flip(traj)
    traj = traj(end:-1:1,:);
end

