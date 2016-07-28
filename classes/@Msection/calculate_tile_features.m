function obj = calculate_tile_features(obj, filter, force, max_features)
% generate features for all tiles in Msection obj
% force = 1: calcualtes features even if the features have already been
% calculated before. Replaces old feature set. default = 0;
% filter can be 'bkgrd1', 'bkgrd2', 'bkgrd3', 'histeq', or 'none' (default), or any
% other filter applicable by tile.get_image

if nargin<2, filter = 'true';end
if nargin<3, force  = 0;end
if nargin<4, max_features = 5000;end
disp('Calculating image features using parfor...');
mL2_tiles = obj.tiles;
%parfor_progress(numel(mL2_tiles));
parfor ix = 1:numel(mL2_tiles)
    %disp(ix);
    if isempty(mL2_tiles(ix).features) || force
         t = get_features(mL2_tiles(ix), filter);
        % reduce feature set (and point set) if too many
        if size(t.features,1)>max_features
            indx = randi(size(t.features,1), max_features,1);
            t.features(indx,:)=[];
            t.validPoints(indx) = [];
        end
        mL2_tiles(ix) = t;
    end
%    parfor_progress;
end
%parfor_progress(0);
disp('Done!');
obj.tiles = mL2_tiles;