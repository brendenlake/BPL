classdef viz_skel
    % visualize the process of extracting the skeleton
    
    methods (Static)
        
        function plot_image(I)
            if UtilImage.check_black_is_true(I)
                I = ~I;
            end
            image(repmat(I,[1 1 3]));
            set(gca,'XTick',[],'YTick',[]);
        end
            
        function plot_junctions(I,J)
           % I : binary image
           % J : junctions
           if UtilImage.check_black_is_true(I)
                I = ~I;
           end
           sz = size(I);
           Z = zeros(sz(1),sz(2),3);
           for i=1:3
               slice = I | J;
               if i>1
                   slice(J) = false;
               end
               Z(:,:,i) = slice;
           end
           image(double(Z));
           set(gca,'XTick',[],'YTick',[]);
        end
        
        function plot_trace(I,T)
           % I : binary image
           % T : graph structure
            if UtilImage.check_black_is_true(I)
                I = ~I;
            end
            sz = size(I);
            hold on
            ns = length(T.S);
            
            image([1 sz(1)],[1 sz(2)],repmat(I,[1 1 3]));
            for i=1:ns
                stk = T.S{i};
                color = rand(3,1);
                plot_traj(stk,color);
            end    
            for i=1:T.n
               plot(T.G(i,2),T.G(i,1),'r.','MarkerSize',18); 
            end
            set(gca,'YDir','reverse','XTick',[],'YTick',[]);
            xlim([1 sz(1)]);
            ylim([1 sz(2)]);
            
        end  
    end
    
end

% plot a stroke trajectory in image space,
% where the x and y dimensions are reversed
function plot_traj(stk,color)
    ystk = stk(:,2);
    stk(:,2) = stk(:,1);
    stk(:,1) = ystk;       
    plot(stk(:,1),stk(:,2),'Color',color,'LineWidth',2);
end