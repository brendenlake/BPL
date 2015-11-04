% Find the smoothest stroke that still fits within the image, and is within a given tolerance of
% error from the original trajectory.
%
% Input
%  stk: [n x 2] trajectory
%  I: binary image
%  bool_viz: 
%  ps: bottom-up parameters
%
function yy = fit_smooth_stk(stk,I,bool_viz,ps)

    slen = size(stk,1); 
    if slen==1 % nothing we can do for points
       yy = stk;
       return
    end
    x = vec(1:slen);
    xx = x;
    count = 1;
    rho = ps.int_rho(count);
    yy = csaps(x,stk',rho,xx);
    yy = yy';
    
    if bool_viz, viz_smoothing(yy,stk,I,rho); end
    
    while ~stk_in_img(yy,I) || ~within_allowed_error(yy,stk,ps)
        count = count + 1;
        if count > numel(ps.int_rho)
           yy = stk;
           break 
        end
        rho = ps.int_rho(count);
        yy = csaps(x,stk',rho,xx);
        yy = yy';
        if bool_viz, viz_smoothing(yy,stk,I,rho); end
    end
    
end

% Are the two trajectories within the allowed error range?
% This is based on the maximum deviation.
function y = within_allowed_error(t1,t2,ps)
    max_dist = compute_dist(t1,t2);
    y = max_dist <= ps.traj_abs_error_lim;
end

% average distance between trajectories
function [max_dist,ave_dist] = compute_dist(t1,t2)
    n1 = size(t1,1);
    n2 = size(t2,1);
    assert(n1==n2)
    dist = zeros(n1,1);
    for i=1:n1
       v1 = t1(i,:);
       v2 = t2(i,:);
       dist(i) = norm(v1(:)-v2(:));
    end
    ave_dist = mean(dist);
    max_dist = max(dist);
end

% Are all the pixels in a stroke "on"?
% This includes the 4 closest pixels to any 
% given point along the line.
%
% Input
%   stk: [n x 2] stroke
%
function bool_in = stk_in_img(stk,I)    
    x = stk(:,1);
    y = stk(:,2);
    
    lb_x = floor(x);
    lb_y = floor(y);
    ub_x = ceil(x);
    ub_y = ceil(y);
    
    lind = [];
    lind = [lind; cut_sub2ind(I,lb_x,lb_y)];
    lind = [lind; cut_sub2ind(I,lb_x,ub_y)];
    lind = [lind; cut_sub2ind(I,ub_x,lb_y)];
    lind = [lind; cut_sub2ind(I,ub_x,ub_y)];
    pix = I(lind);
    bool_in = all(pix);
end

% same as sub2ind, except it 
% removes indices that are too large
% given the image dimensions
function lind = cut_sub2ind(I,x,y)
    sz = size(I);    
    rmvx = x<=0 | x>sz(1);
    rmvy = y<=0 | y>sz(2);
    rmv = rmvx | rmvy;
    x(rmv) = [];
    y(rmv) = [];
    lind = sub2ind(sz,x,y);
end

function viz_smoothing(yy,stk,I,rho)
    S = {{yy}};
    S = space_img_to_motor(S);
    figure(103);
    clf
    plot_motor_to_image(I,S,true,1);
    dist = compute_dist(yy,stk);
    xlabel(['\rho = ',num2str(rho,3)]);
    ylabel(['dist = ',num2str(dist,3)]);
    input('press enter to continue\n','s');
end