%
% Compute a bottom-up character skeleton.
% This algorithm should be deterministic.
%
% Input
%  I: [N x N] binary image (true = black)
%  bool_viz: visualize the results? (default = true)
%
%  Z: [struct] graph structure
%    fields
%    .n: number of nodes
%    .G: [n x 2] node coordinates
%    .E: [n x n boolean] adjacency matrix
%    .EI: [n x n cell], each cell is that edge's index into trajectory list
%         S. It should be a cell array, because there could be two paths
%    .S: [k x 1 cell] edge paths in the image 
%
function U = extract_skeleton(I,bool_viz)
    
    if ~exist('bool_viz','var')
        bool_viz = false;
    end
    
    assert(UtilImage.check_black_is_true(I));
    
    T = make_thin(I); % get thinned image
    J = extract_junctions(T); % get endpoint/junction features of thinned image
    U = trace_graph(T,J,I); % trace paths between features
    B = U.copy();
    U.clean_skeleton;
    
    if bool_viz
       sz = [313 316]; % figure size      
       h = figure;
       pos = get(h,'Position');
       pos(3:4) = sz;
       set(h,'Position',pos);
       
       % visualize original image
       subplot(2,2,1);
       viz_skel.plot_image(I);
       title('Image');
       
       % visualize thinned image
       subplot(2,2,2);
       viz_skel.plot_junctions(T,J); 
       title('Thinned');
       
       % visualize paths in thinned image
       subplot(2,2,3);
       B.plot_skel;
       title('Graph (raw)');
       
       subplot(2,2,4);
       U.plot_skel;
       title('Graph (cleaned)');
       
       set(gcf,'Position',pos);
       pause(.01);
       drawnow
    end
    
end

% Apply thinning algorithm. First it closes holes
% in the image.
%
% Input
%  I: [n x n boolean] raw image.
%    images are binary, where true means "black"
%
% Output
%  T: [n x n boolean] thinned image.
function T = make_thin(I)
    I = bwmorph(I,'fill');
    T = bwmorph(I,'thin',inf);
end