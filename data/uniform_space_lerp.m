% Convert a stroke [x,y] such that it is uniformly sampled in space.
%
% Input
%   stk : [n x 2] stroke
%   dint : target distance between poitns
%
% Output
%   yi : [m x 2] interpolated stroke
%

%%
function stk_yi = uniform_space_lerp(stk,dint)

    %% return if stroke is too short
    n = size(stk,1);
    if n==1
       stk_yi = stk;
       return;
    end
    
    %% compute distance between each point
    dist = zeros(n,1);
    tormv = false(n,1);
    for i=2:n       
        x1 = stk(i,:);
        x2 = stk(i-1,:);
        dist(i) = norm(x1-x2);
        tormv(i) = dist(i) < 1e-4;
    end
    
    %% remove points that are too close
    dist(tormv) = [];
    stk(tormv,:) = [];
    
    %% return if stroke is too short
    n = size(stk,1);
    if n==1
       stk_yi = stk;
       return;
    end
    
    %% cumulative distance
    cumdist = cumsum(dist);
    start_dist = cumdist(1);
    end_dist = cumdist(end);
    x = cumdist(:);
   
    %% 
    nint = round(end_dist/dint);
    nint = max(nint,2);    
    xi = linspace(start_dist,end_dist,nint);
    stk_yi = interp1(x,stk,xi);   
    
    bool_viz = false;
    if bool_viz
       figure
       hold on       
       plot(stk_yi(:,1),stk_yi(:,2),'b.','MarkerSize',18); 
       plot(stk(:,1),stk(:,2),'r.','MarkerSize',6); 
       xlim([0 105]);
       ylim([-105 0]);
    end
    

end