% Extract all the paths between nodes.
%
% Input
%  T: [n x n boolean] thinned image.
%    images are binary, where true means "black"
%  J: [n x n boolean] critical features (junctions
%    and endpoints)
%  I: [n x n boolean] original image
%
% Output
%  U: graph structure
%    fields
%    .n: number of nodes
%    .G: [n x 2] node coordinates
%    .E: [n x n boolean] adjacency matrix
%    .EI: [n x n cell], each cell is that edge's index into trajectory list
%         S. It should be a cell array, because there could be two paths
%    .S: [k x 1 cell] edge paths in the image 
function U = trace_graph(T,J,I)
    
    Z = struct;
    Z.T = T;
    Z.J = J;
    Z.U = false(size(T)); % which pixels have been explained?
    clear T J    
    
    % Get fork points for the junctions
    jindx = find(Z.J);
    [jx jy] = find(Z.J);
    Z.NF = getforks(jx,jy,Z.T);   
    
    % Set up adjacency graph
    S = [];
    n = length(jindx);
    G = [jx jy];
    E = false(n);
    EI = cell(n);
    
    for i=1:n % for each feature
       
       % starting node
       start = [jx(i) jy(i)]; 
       
       % try all branches
       blist = Z.NF{i}; % get all branching options
       nb = size(blist,1);
       for k=1:nb
          
          br = blist(k,:);
          
          if ~Z.U(lind(br,Z)) % If this branch isn't yet explained
             path = [start; br];
             if Z.J(lind(br,Z))
                % criteria so each junction-to-junction is added once
                if lind(br,Z) > lind(start,Z), continue; end
             else
                Z.U(lind(br,Z)) = true;
                tabu = false(size(Z.T));
                tabu(lind(blist,Z)) = true; % tabu list, of all
                    % other branch points from that junction
                tabu(lind(start,Z)) = true;
                [path,Z] = continuepath(path,Z,tabu);    
             end
          
             % Update adjacency matrix
             S = [S; {path}];
             finish = path(end,:);
             snode = i;
             fnode = find(jx==finish(1) & jy==finish(2));
             E(snode,fnode) = true;
             E(fnode,snode) = true;
             eiset1 = [EI{snode,fnode}; length(S)];
             if ~isequal(sort(eiset1(:),1,'ascend'),unique(eiset1(:)))
                error('weird'); 
             end
             eiset2 = [EI{snode,fnode}; length(S)];
             EI{snode,fnode} = eiset1;
             EI{fnode,snode} = eiset2;
             
          end
       end
    end
    
    % Return graph structure
    U = UGraph(G,E,EI,S,I);
end

% Make linear index
%
% Input
%  pts: [n x 2] points to convert (Rows)
function y = lind(pts,Z)
    y = sub2ind(size(Z.T),pts(:,1),pts(:,2));
end

% Get the neighbors for each fork point
function NF = getforks(jx,jy,T)
    n = length(jx);
    NF = cell(n,1);
    for i=1:n
       junc = [jx(i) jy(i)];
       NF{i} = getneighbors(junc,T);  
    end
end

% Continue following the path until you reach a junction
function [path,Z] = continuepath(path,Z,tabu)    
    next = path(end,:);
    while ~Z.J(lind(next,Z))
        next = pathstep(next,Z,tabu);
        Z.U(lind(next,Z)) = true;
        path = [path; next];
        
        % Make sure we can do a loop and return
        % to the start state
        tabu = false(size(tabu));
    end
    
    % Reset any junction to be true
    Z.U(Z.J) = false;
end

% Return all neighbors in the 3 x 3 grid,
% where present nodes are indicated in V
% 
function [neighbors,num] = getneighbors(indx,V)
    i=indx(1);
    j=indx(2);    
    sz = size(V);
    iteri = max(i-1,1):1:min(i+1,sz(1));
    iterj = max(j-1,1):1:min(j+1,sz(2));
    subI = V(iteri,iterj);
    subI(iteri==i,iterj==j) = false;    
    [onx,ony] = find(subI);
    neighbors = [vec(iteri(onx)) vec(iterj(ony))];
    num = size(neighbors,1);
end

% Find the next step in the path
function next = pathstep(curr,Z,tabu)

    [next,num] = getneighbors(curr,Z.T & ~Z.U & ~tabu);
    if num==0
       error('no neighbors to continue path'); 
    end
    
    % if there is more than on possibility
    if num>1
        
       % if possible, pick a joint point
       ind = lind(next,Z);
       isj = Z.J(ind);
       if any(isj) 
          next = next(isj,:);
       end
       
       % pick the closest candidate left
       dists = pdist2(curr,next,'cityblock');
       jindx = argmin(dists);
       next = next(jindx,:);
    end
    
end