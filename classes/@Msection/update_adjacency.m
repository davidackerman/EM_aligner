function obj = update_adjacency(obj)
% update tile adjacency matrix
% Usage: obj = update_adjacency(obj)
%
% Note: 
%       - depends on statistics toolbox because of pdist2
% Author: Khaled Khairy Janelia Research Campus -- 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LARGE = 50000;
if numel(obj.tiles)<LARGE
obj.update_tile_info_switch = -1; % just use the first tile's W and H for all
obj = update_tile_info(obj);
H = obj.tiles(1).H;
W = obj.tiles(1).W;
obj = update_XY(obj);
a = [obj.X(:) obj.Y(:)];

d = pdist2(a,a);        % depends on statistic toolbox  -------- Sosi: not good for large numbers of tiles
dthresh = sqrt(H^2 + W^2) * obj.dthresh_factor;   % diagonal of the tile times factor
obj.A = sparse(triu(d<dthresh,1));
elseif isa(obj.G,'graph')
    obj.A = adjacency(obj.G);
else
    obj.A = [];
    warning('No adjacency matrix was generated');
end

