%
% Assert that we have a valid graph structure M
%
function assert_valid_graph(M)
    
    n = M.n;
    
    % check that we have the correct sizes 
    % of adjacency matrices
    assert(size(M.G,1)==n);
    assert(size(M.E,1)==n);
    assert(size(M.EI,1)==n);
    
    % make sure all of the paths in S
    % have a corresponding edge
    v = tovec(M.EI);
    uq = unique(v);
    assert(isequal(uq(:),(1:length(uq))'));    
    
    % check on the paths throught the graph
    for i=1:size(M.EI,1) % for each edge in the graph
        for j=1:size(M.EI,2)
            
            v = M.EI{i,j}; % get the associated path
            if isempty(v)
               assert(M.E(i,j)==false);
            else
               assert(M.E(i,j)==true);
            end
            
            for k=1:length(v)
               indx = v(k);
               
               % check that the start and end points
               % align with the coordinates stored in G
               mystart =M.S{indx}(1,:);
               myend = M.S{indx}(end,:);               
               assert(aeq(mystart,M.G(i,:)) || aeq(mystart,M.G(j,:)));
               assert(aeq(myend,M.G(i,:)) || aeq(myend,M.G(j,:)));
            end
            
        end
    end
    
    % check that the matrix E and EI are in correspondence
    for i=1:size(M.EI,1) % for each edge in the graph
        for j=1:size(M.EI,2)
            if M.E(i,j)
               assert(~isempty(M.EI(i,j))); 
            end
        end
    end
    
end