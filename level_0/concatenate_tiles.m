function mL3 = concatenate_tiles(mL3_, lambda)
%% % concatenate mL3
if nargin>=2,
    options.lambda = lambda;
else
    options = [];
end
tiles = [];
for pix = 1:numel(mL3_)
    if numel(mL3_(pix).tiles)>2
        tiles = [tiles mL3_(pix).tiles];
    end
end
mL3 = Msection(tiles);
[mL3, A, S] = filter_based_on_tile_area(mL3, options);