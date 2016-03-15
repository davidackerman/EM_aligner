function [L_vec,a] = reduce_to_connected_components(obj)
% obj can be a collection of tiles that are not in one section (i.e. with or without same z-value).
% Outputs an array of Msection objects (also not necessarily limited to one
% z-value) that are contiguous (connected -- in 2 or 3 dimensions) by point-matches in a "graph"
% sense.
% obj is an Msection object that already has a complete pm struct
% The idea is that each connected component can be solved using the matrix
% solver separately. Note: These connected components need to be assembled
% afterwards.
%
%
% Author: Khaled Khairy. Copyright 2016. Janelia Research Campus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% generate a graph object, based on point-matches
obj.G = graph(obj.pm.adj(:,1), obj.pm.adj(:,2), obj.pm.np, {obj.tiles(:).renderer_id});
b = conncomp(obj.G, 'OutputForm', 'vector');  % generate the logical clusters (connected components)
bins = unique(b);
%% 
parfor bix = 1:numel(bins)
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
