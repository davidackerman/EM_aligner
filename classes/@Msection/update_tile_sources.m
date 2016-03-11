function obj = update_tile_sources(obj, owner, project, stack, server)
% updates information source (renderer database) for the tile
% example:
% stack = 'v7_acquire_LC';
% owner='flyTEM';
% project='FAFB00';
% server='http://tem-services.int.janelia.org:8080/render-ws/v1';
% L = update_tile_sources(L, owner, project, stack, server);
%%%%%%%%%%%%%%%%%%%%%
tiles = obj.tiles;
parfor tix = 1:numel(tiles)
    tiles(tix).owner = owner;
    tiles(tix).project = project;
    tiles(tix).stack = stack;
    tiles(tix).server = server;
end
obj.tiles = tiles;