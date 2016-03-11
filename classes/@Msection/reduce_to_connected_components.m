function [L_vec,a] = reduce_to_connected_components(obj)
% This is only for Msection objects that already have a complete pm struct

% generate a graph object, based on point-matches
obj.G = graph(obj.pm.adj(:,1), obj.pm.adj(:,2), obj.pm.np, {obj.tiles(:).renderer_id});
b = conncomp(obj.G, 'OutputForm', 'vector');  % generate the logical clusters (connected components)
bins = unique(b);
%% can be parfor'd ---- sosi --- left for to profile
for bix = 1:numel(bins)
    %disp(['Isolating component: ' num2str(bix)]);
    %     L_vec(bix) = Msection(obj.tiles(b==bins(bix)));
    indx = b==bins(bix);
    if sum(indx)>1
        L_vec(bix) = reduce_to_tile_subset(obj, find(indx));
    else
        L_vec(bix) = Msection(obj.tiles(indx));
    end
    ntiles(bix) = sum(indx);
    
end
[a, c] = sort(ntiles, 'descend');
L_vec = L_vec(c);


% % find the largest cluster
% largest = 1;
% maxbin = 0;

% for ix = 2:numel(bins)
%     if sum((b==ix))>maxbin, maxbin = sum(b==ix);largest = ix;end
% end
% l = reduce_to_tile_subset(obj, find(b==largest));