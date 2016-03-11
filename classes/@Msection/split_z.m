function L_vec_ = split_z(obj)
%%% returns a vector of Msection objects, one for every z value detected within the tiles in obj
unqt = unique([obj.tiles(:).z]);
for lix = 1:numel(unqt)
    z = unqt(lix);
    indx = find([obj.tiles(:).z]==z);
    tiles = obj.tiles(indx);
    L_vec_(lix) = Msection(tiles);
end