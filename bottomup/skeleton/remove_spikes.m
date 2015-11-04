%
% Remove a path between two nodes,
% where one is a singleton (only one edge), and the 
% path is shorter than the ink width (set as a parameter)
%
function T = remove_spikes(T)

    ps = defaultps_bottomup();
    dist_rmv = ps.min_spike_distance; % remove spikes that are shorter than this, 
                  % which is approximately the ink width (like 4)

    degree = sum(T.E);
    degree1 = find(degree==1);
    
    tormv_node = false(T.n,1);
    tormv_edge = false(size(T.S));
    for i=1:length(degree1)
       node1 = degree1(i);
       node2 = find(T.E(node1,:));
       if degree(node2)==1
           continue % don't remove two singletons conjoined 
       end
       edges = T.EI{node1,node2};
       
       if length(edges)==1
           s_edge = T.S{edges};
           dist = norm(s_edge(1,:)-s_edge(end,:));
           if dist < dist_rmv
              tormv_node(node1) = true;
              tormv_edge(edges) = true;
              % fprintf(1,'remove.\n');
           end
           
       end        
    end
    
    T.n = T.n - sum(tormv_node);
    T.G(tormv_node,:) = [];
    T.E(tormv_node,:) = [];
    T.E(:,tormv_node) = [];
    T.EI(tormv_node,:) = [];
    T.EI(:,tormv_node) = [];
    T.S(tormv_edge) = [];
    
    % re-label the edges
    frmv_edge = find(tormv_edge);
    for i=1:numel(T.EI)
       v = T.EI{i};
       for j=1:length(v)
           T.EI{i}(j) = T.EI{i}(j) - sum(frmv_edge < T.EI{i}(j));
       end       
    end
    
    assert_valid_graph(T);
end