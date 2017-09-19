function [Lo, GS] = get_gray_values(Lo, range)
% given Msection object L with pm field as obtained from
% get_section_point_matches, re-samples to produce anew M and
% produces gray values corresponding
% to those points in M
% Note: the convention is to use tiles from "acquire", i.e. post-lens and scale
% Important is that point-matches correspond to the gray-scale image

%% obtain rectangle to sample from
im1 = mat2gray(get_image(Lo.tiles(Lo.pm.adj(1,1))), range);
npoints = 5;
for pix = 1:size(Lo.pm.M,1)
    % determine transformation
    tform = estimateGeometricTransform(Lo.pm.M{pix,1}, Lo.pm.M{pix,2}, 'affine');
    % rectangle determined from first tile
    R1 = round( [min(Lo.pm.M{pix,1}(:,1)) max(Lo.pm.M{pix,1}(:,1))  ...
        min(Lo.pm.M{pix,1}(:,2)) max(Lo.pm.M{pix,1}(:,2))]);
    
    
    x = randi([R1(1), R1(2)],npoints,1);
    y = randi([R1(3), R1(4)],npoints,1);
    X1 = [x(:) y(:)];
    % transform those points to their counterparts in tile 2
    X2 =  transformPointsForward(tform, X1);
    % filter out points outside the image bounds
    indxdel = find(X2(:,1)>size(im1,1));
    indxdel = [indxdel(:); find(X2(:,2)>size(im1,2))];
    indxdel = [indxdel(:); find(X1(:,1)>size(im1,1))];
    indxdel = [indxdel(:); find(X1(:,2)>size(im1,2))];
    X2(indxdel,:) = [];
    X1(indxdel,:) = [];
    Lo.pm.M{pix,1} = [X1(:,1) X1(:,2)];
    Lo.pm.M{pix,2} = [X2(:,1) X2(:,2)];
    %%%% sosi
%     x1 = Lo.pm.M{pix,1};
%     x2 = Lo.pm.M{pix,2};
%     [min(x1) max(x1); min(x2) max(x2)]
    %%%%%%%%%%%%%%%%%%%%%%%%%
end
%% % sosi
    d = 10;
    pix = 1;
x1 = Lo.pm.M{pix,1};
x2 = Lo.pm.M{pix,2};
disp([min(x1) max(x1)]);
disp([min(x2) max(x2)]);
disp(size(im1));
% %
im1 = mat2gray(get_image(Lo.tiles(Lo.pm.adj(pix,1))), range);
im1 = medfilt2(im1, [d d]);
im2 = mat2gray(get_image(Lo.tiles(Lo.pm.adj(pix,2))), range);
im2 = medfilt2(im2, [d d]);

figure;imshow(mat2gray(im1));hold on;plot(Lo.pm.M{pix,1}(:,2),...
    Lo.pm.M{pix,1}(:,1), 'y*');
figure;imshow(mat2gray(im2));hold on;plot(Lo.pm.M{pix,2}(:,2),...
    Lo.pm.M{pix,2}(:,1), 'y*');
figure;showMatchedFeatures(im1, im2, x1, x2, 'montage');
drawnow;
%%

%% determine gray-scale values
d = 10;
GS  = {};
for pix = 1:size(Lo.pm.M,1)
    im1 = mat2gray(get_image(Lo.tiles(Lo.pm.adj(pix,1))), range);
    im1 = medfilt2(im1, [d d]);
    
    indx1 = sub2ind(...
        size(im1), ...
        round(Lo.pm.M{pix,1}(:,1)),...
        round(Lo.pm.M{pix,1}(:,2)));
    
    gs1 = im1(indx1);
    
    im2 = mat2gray(get_image(Lo.tiles(Lo.pm.adj(pix,2))), range);
    im2 = medfilt2(im2, [d d]);
    
    indx2 = sub2ind(size(im2), ...
        round(Lo.pm.M{pix,2}(:,1)),...
        round(Lo.pm.M{pix,2}(:,2)));
    
    gs2 = im2(indx2);
    
    % find and delete entries that point to zero gray scale
    ixdel = find(gs1==0);
    ixdel = [ixdel(:); find(gs2==0)];
    gs1(ixdel) = [];
    gs2(ixdel) = [];
    Lo.pm.M{pix,1}(ixdel,:) = [];
    Lo.pm.M{pix,2}(ixdel,:) = [];
    %
    
    GS{pix,1} = gs1;
    GS{pix,2} = gs2;
    
    
end



























