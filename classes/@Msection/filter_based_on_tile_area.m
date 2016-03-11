function [obj, A, S] = filter_based_on_tile_area(obj, options)
% consider merging with detect_spurious_tiles.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if nargin<2, 
%     areathresh = 1e-2;perimthresh = 1e-2;transthresh = 1e-2;
% else
%     areathresh = options.areathresh;perimthresh = options.perimthresh;
% end
if isfield(options, 'lambda');
lambda = options.lambda;
else
    lambda = 1.0;
end
A = [];
S = [];
trans = [];
obj = update_tile_info(obj);		% we don't need adjacency information in this file

%%% determine polygonal areas
parfor ix = 1:numel(obj.tiles)
    if strcmp(class(obj.tiles(ix).tform), 'affine2d')
        x = obj.tiles(ix).tform.T(3,1);
        y = obj.tiles(ix).tform.T(3,2);
    else
        x = 0;
        y = 0;
    end
    Px = [x; x + obj.tiles(ix).W/obj.map_display_fac; x + obj.tiles(ix).W/obj.map_display_fac; x];
    Py = [y; y; y + obj.tiles(ix).H/obj.map_display_fac; y+obj.tiles(ix).H/obj.map_display_fac];
    
    %%% transform the points
    if strcmp(class(obj.tiles(ix).tform), 'affine2d')
        P = [Px(:) Py(:) [1 1 1 1]']*obj.tiles(ix).tform.T;
    else
        P = transformPointsInverse(obj.tiles(ix).tform,[Px Py]);
    end
    % check polygon area
    A(ix) = polyarea(P(:,1), P(:,2));
    % add polygonperimeter
    s = 0;
    s = s + sqrt((P(1,1)-P(2,1)).^2 + (P(1,2)-P(2,2)).^2);
    s = s + sqrt((P(2,1)-P(3,1)).^2 + (P(2,2)-P(3,2)).^2);
    s = s + sqrt((P(3,1)-P(4,1)).^2 + (P(3,2)-P(4,2)).^2);
    s = s + sqrt((P(1,1)-P(4,1)).^2 + (P(1,2)-P(4,2)).^2);
    S(ix) = s;
    % translation
    %trans(ix,:) =[x y]; 
end

[mu,sig] = normfit(S);

indx = [find(S<(mu-lambda*sig)) find(S>(mu+lambda*sig))];
for ix = 1:numel(indx)
    disp(['Outlier tile found: ' num2str(indx(ix))]);
    obj.tiles(indx(ix)).state = -3;
end

