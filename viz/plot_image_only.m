%
% Plot a motor trajectory on top of the image
%
% Input
%   I [105 x 105] image (double or binary) (large numbers are BLACK)
%   drawings: [nested cell] of strokes in motor space
function plot_image_only(I)
    assert(size(I,1)==105);
    assert(size(I,2)==105);
    newI = 1-I;
    newI(isinf(I)) = 0.5;
    image([1 105],[1 105],repmat(newI,[1 1 3]));
    set(gca,'YDir','reverse','XTick',[],'YTick',[]);
    xlim([1 105]);
    ylim([1 105]);
end