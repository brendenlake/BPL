%
% Removes graph nodes that are within one pixel of eachother
%
% Input
%   T: graph struct that is the output of the tracer
%
% Output
%   T:
%
function T = graph_merge(T)

    % make all merge moves
    tormv = true;
    while tormv
       [T,tormv] = remove_conncomp(T); 
    end

end

% make one removal
function [M,tormv] = remove_conncomp(M)
    
    G = M.G;
    D = squareform(pdist(G));
    EE = D < 1.1; % make an edge between two nodes that are very close
    [S,C] = graphconncomp(sparse(EE),'Directed',false);
    % C indicates which component each node belongs to
    
    tormv= false;
    for i=1:S
       sel = C==i;
       ns = sum(sel);
       if ns > 1
          tormv = true;
          break
       end
    end
    
    % break when necessary
    if ~tormv, return; end
   
    %E_sel_to_sel = M.E(sel,sel);
    E_sel_to_not = M.E(sel,~sel);
    EI_sel_to_sel = M.EI(sel,sel);
    EI_sel_to_not = M.EI(sel,~sel);
    
    % create the new row in adjacency matrix, which 
    % includes all the nodes connected to the selected set/pair
    new_row_E = any(E_sel_to_not,1);
    new_row_EI = cell(size(new_row_E));
    for i=1:size(EI_sel_to_not,2);
        new_row_EI{i} = tovec(EI_sel_to_not(:,i));
    end
 
    % list_ei_sel_to_not = tovec(EI_sel_to_not);
    list_ei_sel_to_sel = tovec(EI_sel_to_sel);
    list_ei_sel_to_sel = unique(list_ei_sel_to_sel);
    
    % see which of the edges, from sel to sel, should be removed.
    % we don't want to remove edges that make a long arc and then return
    nss = length(list_ei_sel_to_sel);
    rmv_edge = false(nss,1);
    for i=1:nss
       indx = list_ei_sel_to_sel(i);
       edge_length = size(M.S{indx},1);
       if edge_length == 2
          rmv_edge(i) = true; 
       end
    end
    keep_ei_sel_to_sel =  list_ei_sel_to_sel(~rmv_edge);
    rmv_ei_sel_to_sel = list_ei_sel_to_sel(rmv_edge);
   
    % remove node positions
    G_sel = M.G(sel,:);
    new_node = mean(G_sel,1);
    Grmv = M.G(sel,:);
    M.G(sel,:) = [];
    M.G = [M.G; new_node];
    M.n = M.n - sum(sel) + 1;
    
    % modify adjacency matrix with the new row
    M.E(sel,:) = [];
    M.E(:,sel) = [];
    M.E(end+1,:) = new_row_E;
    M.E(:,end+1) = [new_row_E,false];
    if ~isempty(keep_ei_sel_to_sel)
        M.E(end,end) = true;
    end
    
    % modify the list of edges with each cell
    M.EI(sel,:) = [];
    M.EI(:,sel) = [];
    M.EI(end+1,:) = new_row_EI;
    M.EI(:,end+1) = [new_row_EI, {keep_ei_sel_to_sel} ];  
    
    % modify the paths
    M.S(rmv_ei_sel_to_sel) = [];
    for i=1:numel(M.EI)
       mycell = M.EI{i};
       
       % update counts in adjacency matrix 
       for j=1:length(mycell)
          el = mycell(j);
          M.EI{i}(j) = el - sum(rmv_ei_sel_to_sel < el);
       end
    end
    
    % update the paths, stored in M.S, such that all 
    % nodes we have now replaced are updated with their
    % new coordinates
    for i=1:length(M.S)
        for j=1:size(Grmv,1);
            node_rmv = Grmv(j,:);
            d = pdist2(node_rmv,M.S{i});
            swap = d < .001;
            nswap = sum(swap);
            M.S{i}(swap,:) = repmat(new_node,[nswap 1]);
        end
    end
    
    assert_valid_graph(M);   
end