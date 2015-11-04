% Convert a stroke [x,y,t] such that it is uniformly sampled
% in time. This is done by linear interpolation in time
%
% Input
%  stk:  [n x 2] for x and y
%  time: [n x 1] for time coordinate
%  tint: [scalar] time interval (in milliseconds)
%
% Output
%  unif_stk: [k x 2] new stroke x and y 
%  unif_t: [k x 1] uniform time interval
%
function [unif_stk,unif_t] = uniform_time_lerp(stk,time,tint)

    mint = min(time);
    maxt = max(time);
    
    unif_t = mint:tint:maxt;
    % make sure we don't leave out the last point
    if unif_t(end) ~= maxt 
       unif_t = [unif_t maxt]; 
    end    
    nt = length(unif_t);
    
    unif_stk = zeros(length(unif_t),2);
    for i=1:nt
        ti = unif_t(i);
        
        diff = time-ti;
        if any(diff==0)
            sel = stk(diff==0,:);            
            unif_stk(i,:) = mean(sel,1);
            continue
        end
        
        % find the point before and after in time
        indx_gt = find(diff>0,1,'first');
        indx_lt = find(diff<0,1,'last');
        
        x_gt = stk(indx_gt,:);
        x_lt = stk(indx_lt,:);
        t_gt = time(indx_gt);
        t_lt = time(indx_lt);        
        
        % Compute the linear interpolation
        frac = (ti - t_lt) ./ (t_gt - t_lt);
        assert(frac<=1);
        unif_stk(i,:) = (1-frac).*x_lt + frac.*x_gt;
    end

end