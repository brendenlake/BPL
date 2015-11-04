%
% Partition a stroke in to sub-strokes based on pauses of the pen.
%
% Input
%  unif_stk: [n x 2] stroke, assuming uniform time sampling
%  dthresh: [scalar] if this much distance (norm) is not covered
%         at each time point, then it's a pause
%  max_sequence: [scalar] maximum length of a stop sequence, before it is
%         called it's own stroke
%
% Output
%  substrokes: [ns x 1] cell array of sub-strokes
%  unif_stk: stroke with pause sequences shortened to a single point
%  breaks: where the pauses occured
%
function [substrokes,unif_stk,breaks] = partition_strokes(unif_stk,dthresh,max_sequence)

    % Compute derivatives
    n = size(unif_stk,1);
    dxdt = get_deriv(unif_stk,1:n);    
    
    % Special case
    if n==1
       substrokes = {unif_stk};
       breaks = true;
       return
    end
    
    % Compute norm of derivatives
    norm_dxdt = zeros(n,1);
    for i=1:n
        norm_dxdt(i) = norm(dxdt(i,:)); 
    end
    
    %% 
    bool_viz = false;
    if bool_viz
       figure(99);
       clf
       title('histogram of distances');
       hist(norm_dxdt,round(n/2));       
    end    
    
    %% Compute the candidate stop points
    stop_pt = norm_dxdt < dthresh;
    for i=2:n
       if stop_pt(i)
           stop_pt(i-1) = true;
       end
    end
    stop_pt(1) = true;
    stop_pt(end) = true;

    %% 
    if bool_viz
       figure(100);
       clf
       hold on
       plot(unif_stk(:,1),unif_stk(:,2),'b.','MarkerSize',12);
       plot(unif_stk(stop_pt,1),unif_stk(stop_pt,2),'g.','MarkerSize',12);
    end 
    
    % Partition the stop points into stop sequences.
    % Here, non-stops are denoted as zeros, the first stop
    % is a sequence of 1s, second is a sequence of twos, etc.
    % Until the pen is moving fast enough again
    stop_sequence = zeros(n,1);
    stop_count = 1;
    for i=1:n        
        if stop_pt(i) % current point is a stop, it's the same stop
            stop_sequence(i) = stop_count;
        elseif stop_pt(i-1) && stop_pt(i+1) % points surround it are a stop... its the same stop
            stop_sequence(i) = stop_count;
        elseif stop_pt(i-1)
            stop_count = stop_count + 1; % we just finishsed a stop
        end
    end
    
    % Special case where the entire stroke is a stop sequence
    if stop_count == 1
        stop_sequence = zeros(n,1);
        stop_sequence(1) = 1;
        stop_sequence(end) = 2;
        stop_count = 2;
    end    
    
    % Make sure the stop sequences aren't too long. If they are,
    % we place a sub-stroke break at the beginning and end.
    i = 1;
    while i <= stop_count
        sel = find(stop_sequence==i);
        nsel = length(sel);
        if nsel>max_sequence
            stop_sequence(sel(2:end)) = 0;
            stop_sequence(stop_sequence>i) = stop_sequence(stop_sequence>i) + 1;
            stop_sequence(sel(end)) = i+1;
            stop_count = stop_count + 1;
        end
        i = i + 1;
    end    
    
    % Breaks are the average of the stop sequences
    mybreaks = zeros(n,1);
    for i=1:stop_count
        sel = stop_sequence==i; % select the stop sequence
        
        if i==1 % beginning of stroke
            mybreaks(i) = find(sel,1,'first');
        elseif i==stop_count % end of stroke
            mybreaks(i) = find(sel,1,'last');
        else % all other positions
            mybreaks(i) = round(mean(find(sel))); % find the mean element
        end
            
        % Set the mean element to the mean of the sequence
        unif_stk(mybreaks(i),:) = mean(unif_stk(sel,:),1);
        
        % mark to keep
        stop_sequence(mybreaks(i)) = -1;
    end
    
    % Remove all other stop sequence elements, 
    % except for the marked mean
    unif_stk(stop_sequence>0,:) = [];
    stop_sequence(stop_sequence>0) = [];
    breaks = stop_sequence<0;
    
    % Convert to cell array
    fbreaks = find(breaks);
    nb = length(fbreaks);
    ns = max(1,nb-1);
    substrokes = cell(ns,1);
    if nb==1 % if this stroke was just a single stop sequence
        assert(size(unif_stk,1)==1);
        substrokes{1} = unif_stk;
    else
        for s=1:ns
            substrokes{s} = unif_stk(fbreaks(s):fbreaks(s+1),:);
        end
    end
    
    new_start = substrokes{1}(1,:);
    new_end = substrokes{end}(end,:);
    assert( aeq(new_start,unif_stk(1,:)) );
    assert( aeq(new_end,unif_stk(end,:)) );
end

% Compute dx/dt partial derivatives
% for each column (variable) of X
% 
% Input
%  X: [T x dim]
%  t: [T x 1] time points
% 
% Output
%  dxdt: [T-1 x 1] derivatives
function dxdt = get_deriv(X,t)

    [T,dim] = size(X);
    assert(isvector(t));
    assert(numel(t)==T);
    dxdt = zeros(size(X));
    
    for i=2:T
        prev = X(i-1,:);
        next = X(i,:);
        dt = t(i)-t(i-1);
        dxdt(i,:) = (next-prev)./dt;
    end
    
end