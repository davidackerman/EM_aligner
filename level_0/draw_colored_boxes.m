function  draw_colored_boxes(boxes, data, data_bounds, title_str)
% generate figure with colored boxes for diagnostics
% used by gen_diagnostics.m
if nargin<4, title_str = 'untitled figure';end
figure;
c = mat2gray(data, data_bounds);

for tix = 1:numel(boxes)
    P = boxes{tix}{1};   % patch for this tile
    patch( P(:,1), P(:,2), [c(tix)],'FaceColor', 'flat',   'EdgeColor', 'none' , 'Facealpha', 0.4);
end
title(title_str);
daspect([1 1 1]);axis ij; h = colorbar; set(h, 'ylim', data_bounds);
drawnow;