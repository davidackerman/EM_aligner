function [Lo, GS] = get_gray_values(Lo, range)
% given Msection object L with pm field as obtained from
% get_section_point_matches, re-samples to produce anew M and 
% produces gray values corresponding 
% to those points in M
% Note: the convention is to use tiles from "acquire", i.e. post-lens and scale
% Important is that point-matches correspond to the gray-scale image

%% obtain rectangle to sample from
npoints = 50;
for pix = 1:size(Lo.pm.M,1)
    % determine transformation
    tform = estimateGeometricTransform(Lo.pm.M{pix,1}, Lo.pm.M{pix,2}, 'affine');
% rectangle determined from first tile
R1 = round( [min(Lo.pm.M{pix,1}(:,2)) max(Lo.pm.M{pix,1}(:,2))  ...
             min(Lo.pm.M{pix,1}(:,1)) max(Lo.pm.M{pix,1}(:,1))]);

         
 x = randi([R1(1), R1(2)],npoints,1);
 y = randi([R1(3), R1(4)],npoints,1);
 X1 = [x(:) y(:)];
 % transform those points to their counterparts in tile 2
 X2 =  transformPointsForward(tform, X1);
 Lo.pm.M{pix,1} = [X1(:,1) X1(:,2)];
 Lo.pm.M{pix,2} = [X2(:,1) X2(:,2)];
end


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
    ixdel = [ixdel(:) find(gs2==0)];
    gs1(ixdel) = [];
    gs2(ixdel) = [];
    Lo.pm.M{pix,1}(ixdel,:) = [];
    Lo.pm.M{pix,2}(ixdel,:) = [];
    %
    
    GS{pix,1} = gs1;
    GS{pix,2} = gs2;
    
    %%% sosi
%     figure;imshow(mat2gray(im1));hold on;plot(Lo.pm.M{pix,1}(:,2),...
%         Lo.pm.M{pix,1}(:,1), 'y*');
%     figure;imshow(mat2gray(im2));hold on;plot(Lo.pm.M{pix,2}(:,2),...
%         Lo.pm.M{pix,2}(:,1), 'y*');
%     showMatchedFeatures(im1, im2, Lo.pm.M{pix,1}, Lo.pm.M{pix,2}, 'montage');
%     drawnow;
end



























