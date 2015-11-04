%
% Convert a random walk to a canonical
% format, so they can be compared.
%
% Input
%  S: stroke array in image space
%
% Each stroke is arranged so it goes "top-to-bottom"
%
% Then the strokes are ordered such they progress
% from top-to-bottom, in terms of starting location
%
function S = walk_to_canonical(S)

    S = space_img_to_motor(S);

    % choose canonical directions
    ns = length(S);
    for i=1:ns
       stk = S{i};
       if distance_TL(stk(end,:)) < distance_TL(stk(1,:))
           S{i} = flip_stroke(S{i});
       elseif is_circle(stk) ... % special case of circular strokes
           && distance_TL(stk(end-1,:)) < distance_TL(stk(2,:))
           S{i} = flip_stroke(S{i});
       end
    end
    
    % compute distance of starting location to corner,
    % or ending location to the corner
    mystart = zeros(ns,2);
    myend = zeros(ns,2);
    for i=1:ns
       mystart(i,:) = S{i}(1,:);
       myend(i,:) = S{i}(end,:);
    end
    d_start = distance_TL(mystart);
    %[~,end_indx] = sort(myend(:,2),1,'ascend');
    
    % break ties by their end location
    for i=1:numel(d_start)
       sel = d_start == d_start(i);
       if sum(sel) > 1 % if we have a tie
            [~,end_indx] = sort(myend(sel,2),1,'ascend'); % sort them
            end_indx = end_indx-1;
            d_start(sel) = d_start(sel) + end_indx/10;
       end
    end
    
    % Sort based on starting location (with ties broken)
    [~,indx] = sort(d_start,1,'ascend');
    S = S(indx);
    S = space_motor_to_img(S);
    S = S(:);
end

% flip a trajectory
function stk = flip_stroke(stk)
    stk = stk(end:-1:1,:);
end

% does the stroke have the same start and end?
function y = is_circle(stk)
    y = size(stk,1) >= 4 && isequal(stk(1,:),stk(end,:));
end

%
% Compute the distance from the top-left
% of each image to the current point
% 
% Only works if everything is in motor space
%
function d = distance_TL(list_pts)
    [n,dim] = size(list_pts);
    assert(dim==2);
    d = zeros(n,1);
    for i=1:n
       d(i) = norm(list_pts(i,:)); 
    end
end

