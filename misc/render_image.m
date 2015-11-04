% 
% Render a gray-scale probability map image, given the sub-stroke splines
%
% Input
%  cell_traj: [nested cell] b-spline coordinates of each sub-stroke
%  epsilon: probability of a pixel flip
%  blur_sigma: bandwidth of blur noise
%  PM: parameters from defaultps
%
% Output
%  prob_on: probability of each pixel being on
%
function [prob_on,ink_off_page] = render_image(cell_traj,epsilon,blur_sigma,PM)
    
    % convert to image space
    traj_img = space_motor_to_img(cell_traj);
    
    % draw the trajectories on the image  
    template = zeros(PM.imsize);
    nsub = length(traj_img);
    ink_off_page = false;
    for i=1:nsub
       
       % Ink model parameters
       ink = PM.ink_pp; % amount of ink per point
       max_dist = PM.ink_max_dist; % distance to which you get full ink
        
       % check boundaries
       myt = traj_img{i};
       out = check_bounds(myt,PM.imsize);
       if any(out), ink_off_page = true; end
       myt(out,:) = [];
       if isempty(myt), continue; end
       
       % compute distance between each trajectory point
       % and the next one
       if size(myt,1)==1
           myink = ink;
       else
           dist = pair_dist(myt);
           dist(dist>max_dist) = max_dist;
           dist = [dist(1); dist];
           myink = ink .* dist ./ max_dist;
       end
       
       % Make sure we have the minimum amount of ink,
       % if a particular trajectory is very small
       sumink = sum(myink);
       if aeq(sumink,0)
          nink = numel(myink); 
          myink = ones(size(myink)) .* ink ./ nink;
       elseif sumink < ink
          myink = (ink ./ sumink) .* myink; 
       end
       assert(sum(myink)>ink-1e-4);       
       
       % share ink with the neighboring 4 pixels
       x = myt(:,1);
       y = myt(:,2);
       xfloor = floor(x);
       yfloor = floor(y);
       xceil = ceil(x);
       yceil = ceil(y);
       x_c_ratio = x - xfloor;
       y_c_ratio = y - yfloor;
       x_f_ratio = 1 - x_c_ratio;
       y_f_ratio = 1 - y_c_ratio;       
       lin_ff = sub2ind(PM.imsize,xfloor,yfloor);
       lin_cf = sub2ind(PM.imsize,xceil,yfloor);
       lin_fc = sub2ind(PM.imsize,xfloor,yceil);
       lin_cc = sub2ind(PM.imsize,xceil,yceil);
       
       % paint the image
       template = seqadd(template,lin_ff,myink .* x_f_ratio .* y_f_ratio);
       template = seqadd(template,lin_cf,myink .* x_c_ratio .* y_f_ratio);
       template = seqadd(template,lin_fc,myink .* x_f_ratio .* y_c_ratio);
       template = seqadd(template,lin_cc,myink .* x_c_ratio .* y_c_ratio);
    end
   
    % filter the image to get the desired brush-stroke size
    a = PM.ink_a;
    b = PM.ink_b;
    H_broaden = b .* [a/12 a/6 a/12; a/6 1-a a/6; a/12 a/6 a/12];
    widen = template;
    for i=1:PM.ink_ncon
        widen = imfilter(widen,H_broaden,'conv');
    end
    
    % threshold again
    widen(widen>1) = 1;   
    
    % filter the image to get Gaussian
    % noise around the area with ink
    pblur = widen;
    if blur_sigma > 0
        fsize = 11; 
        H_gaussian = fspecial('gaussian',fsize,blur_sigma);
        pblur = imfilter(pblur,H_gaussian,'conv');
        pblur = imfilter(pblur,H_gaussian,'conv');
    end

    % final truncation
    pblur(pblur > 1) = 1;
    pblur(pblur < 0) = 0;
    
    % probability of each pixel being on
    prob_on = (1-epsilon) .* pblur + epsilon .* (1-pblur);    
end

% H_broaden = b * (1+a) .* [a/12 a/6 a/12; a/6 1-a a/6; a/12 a/6 a/12]; % Hinton's original filter

%
% Indices may repeat..
%
% x: [n x 1] data
% lind [k x 1] linear indices in x
% inkval [k x 1] to be added to data, at index lind
%
function x = seqadd(x,lind,inkval)
    for i=1:numel(lind)
        x(lind(i)) = x(lind(i)) + inkval(i);
    end
    
    %
    % These two things are not the same!
	% x(lind) = x(lind) + inkval;
end

function out = check_bounds(myt,imsize)
    xt = myt(:,1);
    yt = myt(:,2);
    out = floor(xt) < 1 | ceil(xt) > imsize(1) | floor(yt) < 1 | ceil(yt) > imsize(2);
end

function z = pair_dist(D)
    x1 = D(1:end-1,:);
    x2 = D(2:end,:);
    z = sqrt(sum((x1-x2).^2,2));
end