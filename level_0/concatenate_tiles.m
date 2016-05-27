function mL3 = concatenate_tiles(mL3_, lambda)
%% % concatenate mL3
if nargin>=2,
    options.lambda = lambda;
else
    options = [];
end

% how many tiles total do we have
tilecount = zeros(numel(mL3_),1);
for pix = 1:numel(mL3_)
    tilecount(pix) = numel(mL3_(pix).tiles);
end
tilecount(tilecount<2)=0;
ntiles = sum(tilecount);
% preallocate and fill tiles array
tiles(ntiles) = tile;
cnt = 1;
for pix = 1:numel(mL3_)
    if tilecount(pix)>=2
    tiles(cnt:cnt+tilecount(pix)-1) = mL3_(pix).tiles;
    end
    cnt = cnt + tilecount(pix);
end
mL3 = Msection(tiles);
[mL3, A, S] = filter_based_on_tile_area(mL3, options);