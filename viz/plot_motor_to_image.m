%
% Plot a motor trajectory on top of the image
%
% Input
%   I [105 x 105] image (double or binary) (large numbers are BLACK)
%   drawings: [nested cell] of strokes in motor space
%   size_start : [scalar] size of start marker
%   lw : line widht of stroke
function plot_motor_to_image(I,drawing,size_start,lw)

    if ~exist('size_start','var') || isempty(size_start) 
       size_start = 18;
    end
    bool_start = size_start > 0;
    if ~exist('lw','var') || isempty(lw)
       lw = 4; % line width 
    end

    assert(size(I,1)==105);
    assert(size(I,2)==105);
    drawing = space_motor_to_img(drawing);
    hold on
    plot_image_only(I);
    
    % Draw each stroke
    ns = length(drawing);
    for i=1:ns
        stroke = drawing{i};
        nsub = length(stroke);
        color = get_color(i);
        for b=1:nsub
            plot_traj(drawing{i}{b},color,lw);
        end
    end
    
    % Plot sub-stroke breaks
    if bool_start
        for i=1:ns
            stroke = drawing{i};
            nsub = length(stroke);
            color = get_color(i);
            for b=1:nsub
                plot_break(drawing{i}{b},size_start);
            end
        end
    end
    
    % Plot starting locations and stroke order markers
    if bool_start
        for i=1:ns
            stroke = drawing{i};
            plot_start_loc(stroke{1}(1,:),i,size_start);
        end
    end
      
    set(gca,'YDir','reverse','XTick',[],'YTick',[]);
    xlim([1 105]);
    ylim([1 105]);
end

function plot_traj(stk,color,lw)
    ystk = stk(:,2);
    stk(:,2) = stk(:,1);
    stk(:,1) = ystk;       
    plot(stk(:,1),stk(:,2),'Color',color,'LineWidth',lw);
end

function plot_break(stk,size_start)
    ystk = stk(:,2);
    stk(:,2) = stk(:,1);
    stk(:,1) = ystk;       
    plot(stk(end,1),stk(end,2),'MarkerEdgeColor','k','Marker','o','MarkerFaceColor','w','MarkerSize',size_start/2) 
end

% Plot the starting location of a stroke
%
% Input
%  start: [1 x 2]
%  num: number that denotes stroke order
function plot_start_loc(start,num,sz_start)
    plot(start(2),start(1),'o','MarkerEdgeColor','k','MarkerFaceColor','w',...
        'MarkerSize',sz_start);
    text(start(2),start(1),num2str(num),...
        'BackgroundColor','none','EdgeColor','none',...
        'FontSize',sz_start,'FontWeight','normal','HorizontalAlignment','center');
end

% Color map for the stroke of index k
function out = get_color(k)
    scol = {'r',[0,0.8,0],'b','m',[0,0.8,0.8]};
    ncol = length(scol);
    if k <=ncol
       out = scol{k}; 
    else
       out = scol{end}; 
    end
end