classdef UGraph < BetterHandle
    %
    % UGraph : Undirected graph for an image skeleton
    %    
    properties
        G  % [n x 2] node coordinates
        E  % [n x n boolean] adjacency matrix
        EI % [n x n cell], each cell is that edge's index into trajectory list
           %  S. It should be a cell array, because there could be two paths
        S  % [k x 1 cell] edge paths in the image
        I  % [nimg x nimg boolean] original image
    end
    
    properties (Dependent = true)
        n % number of nodes in a graph
        imsize % image size
        link_ei_to_ni % [k x 2] for each edge (row), the source and dest. node,
                      % given the direction it is listed in S
    end
    
    properties (SetAccess = private)
       list_mask % circular masks for merging critical points
    end
    
    methods
        
        % constructor for making an undirected graph
        function this = UGraph(G,E,EI,S,I)
           this.G = G;
           this.E = E;
           this.EI = EI;
           this.S = S; 
           this.I = I;          
           
           load('circle_masks','list_mask');
           this.list_mask = list_mask;
        end
        
        function y = get.n(this)
           y = size(this.G,1);            
        end
        
        function sz = get.imsize(this)
            sz = size(this.I);            
        end
        
        % get what nodes the edges are connceted to
        function Y = get.link_ei_to_ni(this)
           k = numel(this.S);
           Y = zeros(k,2);
           for i=1:k
               source_pt = this.S{i}(1,:);
               dest_pt = this.S{i}(end,:);
               source_ni = this.map_pts_to_ni(source_pt);
               dest_ni = this.map_pts_to_ni(dest_pt);
               Y(i,1) = source_ni;
               Y(i,2) = dest_ni;
           end
        end
        
        % Plot the skeleton on-top of the image
        function plot_skel(this)
           this.plot_on_image(~this.I);
        end
        
        % Plot the maximum circle analysis, with the
        % original image and the skeleton
        function plot_circle(this)            
            I_cluster = this.max_circle_mask;
            I_mask = double(this.I);
            I_mask(I_cluster) = 0.5;
            this.plot_on_image(1-I_mask);
        end
        
        % Merge groups of graph nodes connected
        % by the maximum circle criterion
        function clean_skeleton(this)
           I_mask = this.max_circle_mask();
           
           % cluster nodes
           C = this.cluster_G(I_mask);
           nc = numel(unique(C));
           
           % until we no longer have things to remove
           while has_rmv(C)
              
              % for each cluster index, look for 
              % things to merge
              for i=1:nc
                  sel = C == i;
                  if sum(sel)>1
                     this.merge(sel,I_mask);
                     
                     % cluster nodes
                     C = this.cluster_G(I_mask);
                     [~,nc] = has_rmv(C);
                     break
                  end                  
              end
               
           end
           
           % do we still have groups to remove?
           %
           % max_c: maximum cluster index
           function [yyyy,max_c] = has_rmv(myC)
               myC(myC==0) = [];
               yyyy = numel(myC) ~= numel(unique(myC));
               max_c = max(myC);
           end           
           
        end
        
        % Image based on maximum circles surrounding
        % each critical point
        function I_mask = max_circle_mask(this)
            I_mask = false(this.imsize);
            for indx=1:this.n
                I_mask = this.max_circle(indx,I_mask);
            end
        end
        
        
        % merge a collection of nodes
        % 
        % Input
        %   sel: [n x 1 logical] or [k x 1 index]
        %
        function merge(this,sel,I_mask)
            
            % process input
            if islogical(sel)
                assert(numel(sel)==this.n);
            else
                choose = false(this.n,1);
                choose(sel) = true;
                sel = choose;
            end
            
            % partition
            E_sel_to_not = this.E(sel,~sel);
            EI_sel_to_sel = this.EI(sel,sel);
            EI_sel_to_not = this.EI(sel,~sel);
            
            % create the new row in adjacency matrix, which 
            % includes all the nodes connected to the selected set/pair
            new_row_E = any(E_sel_to_not,1);
            new_row_EI = cell(size(new_row_E));
            for i=1:size(EI_sel_to_not,2);
                new_row_EI{i} = tovec(EI_sel_to_not(:,i));
            end
            list_ei_sel_to_sel = tovec(EI_sel_to_sel);
            list_ei_sel_to_sel = unique(list_ei_sel_to_sel);
            
            % see which of the edges, from sel to sel, should be removed.
            % we don't want to remove edges that make a long arc and then return
            nss = length(list_ei_sel_to_sel);
            rmv_edge = false(nss,1);
            for i=1:nss
               indx = list_ei_sel_to_sel(i);       
               
               % remove edge if all the pixels are in the circle criterion mask
               pix = this.S{indx};
               sz = size(I_mask);
               ind_px = sub2ind(sz,pix(:,1),pix(:,2));
               if all(I_mask(ind_px));
                    rmv_edge(i) = true; 
               end
               
            end
            keep_ei_sel_to_sel =  list_ei_sel_to_sel(~rmv_edge);
            rmv_ei_sel_to_sel = list_ei_sel_to_sel(rmv_edge);

            % remove node positions
            G_sel = this.G(sel,:);
            new_node = mean(G_sel,1);
            Grmv = this.G(sel,:);
            this.G(sel,:) = [];
            this.G = [this.G; new_node];

            % modify adjacency matrix with the new row
            this.E(sel,:) = [];
            this.E(:,sel) = [];
            this.E(end+1,:) = new_row_E;
            this.E(:,end+1) = [new_row_E,false];
            if ~isempty(keep_ei_sel_to_sel)
                this.E(end,end) = true;
            end
            
            % modify the list of edges with each cell
            this.EI(sel,:) = [];
            this.EI(:,sel) = [];
            this.EI(end+1,:) = new_row_EI;
            this.EI(:,end+1) = [new_row_EI, {keep_ei_sel_to_sel} ];  

            % modify the paths
            this.S(rmv_ei_sel_to_sel) = [];
            for i=1:numel(this.EI)
               mycell = this.EI{i};

               % update counts in adjacency matrix 
               for j=1:length(mycell)
                  el = mycell(j);
                  this.EI{i}(j) = el - sum(rmv_ei_sel_to_sel < el);
               end
            end

            % update the paths, stored in this.S, such that all 
            % nodes we have now replaced are updated with their
            % new coordinates
            for i=1:length(this.S)
                for j=1:size(Grmv,1);
                    node_rmv = Grmv(j,:);
                    d = pdist2(node_rmv,this.S{i});
                    swap = d < .001;
                    nswap = sum(swap);
                    this.S{i}(swap,:) = repmat(new_node,[nswap 1]);
                end
            end
            
            % Make sure all the trajectories are roughly uniform distance
            for i=1:length(this.S)
                this.S{i} = expand_unif_interp(this.S{i},1);
            end
            
            this.assert_valid_graph();
        end
        
        % Check to make sure this is a valid graph structure
        % Throw assert if it is not
        function assert_valid_graph(this)
            
            n = this.n;
    
            % check that we have the correct sizes 
            % of adjacency matrices
            assert(size(this.G,1)==n);
            assert(size(this.E,1)==n);
            assert(size(this.EI,1)==n);

            % make sure all of the paths in S
            % have a corresponding edge
            v = tovec(this.EI);
            uq = unique(v);
            assert(isequal(uq(:),(1:length(uq))'));    

            % check on the paths throught the graph
            for i=1:size(this.EI,1) % for each edge in the graph
                for j=1:size(this.EI,2)

                    v = this.EI{i,j}; % get the associated path
                    if isempty(v)
                       assert(this.E(i,j)==false);
                    else
                       assert(this.E(i,j)==true);
                    end

                    for k=1:length(v)
                       indx = v(k);

                       % check that the start and end points
                       % align with the coordinates stored in G
                       mystart = this.S{indx}(1,:);
                       myend = this.S{indx}(end,:);               
                       assert(aeq(mystart,this.G(i,:)) || aeq(mystart,this.G(j,:)));
                       assert(aeq(myend,this.G(i,:)) || aeq(myend,this.G(j,:)));
                    end

                end
            end

            % check that the matrix E and EI are in correspondence
            for i=1:size(this.EI,1) % for each edge in the graph
                for j=1:size(this.EI,2)
                    if this.E(i,j)
                       assert(~isempty(this.EI(i,j))); 
                    end
                end
            end
                   
        end
        
        %
        % Map list of nodes (as points) to their indices
        %
        % Input
        %  vpts: [n x 2] points
        %
        % Output
        %  vni: [n x 1] indices
        function vni = map_pts_to_ni(this,vpts)
            [n,dim] = size(vpts);
            assert(dim==2);
            vni = zeros(n,1);
            for i=1:n
                pt = vpts(i,:);
                dist = pdist2(pt,this.G);
                assert(isvector(dist));
                dist = dist(:);
                [minval,minindx] = min(dist);
                assert(aeq(minval,0));
                vni(i) = minindx;
            end   
        end
        
       % Given the current point, where can we go?
        %
        % Input: pt: [1x2] current point we want to analyze
        %
        % Output
        %  vei: [m x 1] list of edges we can go to
        %  vei_flip: [m x 1 boolean] do we need to flip that edge
        %    to align it to the current point
        function [vei,vei_flip] = get_branches(this,pt)
            assert(numel(pt)==2);
            ni = this.map_pts_to_ni(pt);
            list_ei = this.EI(ni,:);
            vei = concat_cell(list_ei);
            vei = vei(:);
            nn = numel(vei);
            list_traj = this.S(vei);
            vei_flip = false(size(vei));
            for i=1:nn
               traj = list_traj{i};
               if isequal(pt,traj(1,:))
                   vei_flip(i) = false;
               elseif isequal(pt,traj(end,:))
                   vei_flip(i) = true;
               else
                   error('cannot map point to its paths'); 
               end
            end            
        end              
        
    end
    
    methods (Access = private)
         
        
        % find the maximum circle in the image surrounding a point
        function I_cluster = max_circle(this,indx,I_cluster)
            pt = this.G(indx,:);
            pt = round(pt);
            nmask = numel(this.list_mask);
            bool_fit = false(nmask,1);
            for i=1:nmask
               x_mask = this.list_mask{i};               
               x_mask(:,1) = x_mask(:,1) + pt(1);
               x_mask(:,2) = x_mask(:,2) + pt(2);               
               x_mask = this.cut_image_plane(x_mask);
               xlind = sub2ind(this.imsize,x_mask(:,1),x_mask(:,2));      
               bool_fit(i) = all(this.I(xlind));
               if bool_fit(i)
                   I_cluster(xlind) = true; 
               else
                   return
               end
            end
        end 
        
        % remove pixels that are out of the image plane
        function x = cut_image_plane(this,x)
            sz = this.imsize;
            xx = x(:,1);
            xy = x(:,2);
            rmvx = xx<=0 | xx>sz(1);
            rmvy = xy<=0 | xy>sz(2);
            rmv = rmvx | rmvy;
            x(rmv,:) = [];
        end
        
        %
        %  Plot the graph skeleton ontop of an image
        %
        function plot_on_image(this,I)            
            sz = size(I);
            if size(I,3) == 1
                I = repmat(I,[1 1 3]);
            end            
            hold on
            ns = length(this.S);
            image([1 sz(1)],[1 sz(2)],I);
            for i=1:ns
                stk = this.S{i};
                %color = rand(3,1);
                plot_traj(stk,'g');
            end    
            for i=1:this.n
               plot(this.G(i,2),this.G(i,1),'r.','MarkerSize',7); 
            end
            set(gca,'YDir','reverse','XTick',[],'YTick',[]);
            xlim([1 sz(1)]);
            ylim([1 sz(2)]);
            
        end
        
        % Cluster the vertices in G
        % into clusters based on the maximum circle criterion
        function C = cluster_G(this,I_mask)
            assert(UtilImage.check_black_is_true(I_mask));
            L = bwlabel(I_mask,4);
            GG = round(this.G);
            lind_G = sub2ind(this.imsize,GG(:,1),GG(:,2));
            C = L(lind_G);            
        end
        
        
    end
    
end

% flatten a cell array into a regular array
% by vert cat each cell
function v = concat_cell(vcell)
    assert(iscell(vcell));
    v = [];
    for i=1:numel(vcell)
        v = [v; vcell{i}];
    end
end

% plot a stroke trajectory in image space,
% where the x and y dimensions are reversed
function plot_traj(stk,color)
    ystk = stk(:,2);
    stk(:,2) = stk(:,1);
    stk(:,1) = ystk;       
    plot(stk(:,1),stk(:,2),'Color',color,'LineWidth',1);
end