% Simple functino for visualizing the trajectory
% of a motor program M
function simpleViz(M)

    cell_traj = M.motor_warped;
    ns = numel(cell_traj);
    
    hold on
    for sid=1:ns
        stroke = cell_traj{sid};
        nsub = length(stroke);        
        for b=1:nsub
            plot_traj(cell_traj{sid}{b});
        end        
    end
    
    xlim([0 105]);
    ylim([-105 0]);
    box on
end

function plot_traj(T)
    plot(T(:,1),T(:,2),'k','LineWidth',2);
end