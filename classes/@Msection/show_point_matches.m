function show_point_matches(obj, flag)

% if flag
% M = obj.pm.M;
% adj = obj.pm.adj;
% using_SURF = 0;
% for pix = 1:size(M,1) % loop over point matches
%     %%%%%transform points for the first of the two tiles
%     if using_SURF,
%         pm = M{pix,1}(:).Location;
%     else   pm = M{pix,1};
%     end
%     T = obj.tiles(adj(pix,1)).tform.T;
%     pmt = pm*T(1:2,1:2);
%     if using_SURF,
%         M{pix,1}(:).Location = pmt;
%     else
%         M{pix,1}(:) = pmt;
%     end
%     %%%%%%%%%%%transform points for the second of the two tiles
%     if using_SURF
%         pm = M{pix,2}(:).Location;
%     else
%         pm = M{pix,2};
%     end
%     T = obj.tiles(adj(pix,2)).tform.T;
%     pmt = pm*T(1:2,1:2);
%     if using_SURF
%         M{pix,2}(:).Location = pmt;
%     else
%         M{pix,2} = pmt;
%     end
% end
% obj.pm.M = M;
% obj.pm.adj = adj;
% end


if numel(obj.tiles)==2
   fac = 0.4;
   maxM = 3; 
   m = randi(size(obj.pm.M{1,1},1), maxM,1);
   im1 = imresize(get_image(obj.tiles(1)), fac);
   im2 = imresize(get_image(obj.tiles(2)), fac);
   showMatchedFeatures(im1, im2,obj.pm.M{1,1}(m,:) * fac,obj.pm.M{1,2}(m,:) * fac, 'Method', 'Montage');
    
    
else
show_map(obj);
hold on;
for ix = 1:size(obj.pm.adj,1)
    m1ix = obj.pm.adj(ix,1);
    m2ix = obj.pm.adj(ix,2);
    m1 = obj.pm.M{ix,1};
    m2 = obj.pm.M{ix,2};
    
    
%     x1 = obj.tiles(obj.pm.adj(ix,1)).tform.T([3]);
%     y1 = obj.tiles(obj.pm.adj(ix,1)).tform.T([6]);
%     
%     x2 = obj.tiles(obj.pm.adj(ix,2)).tform.T([3]);
%     y2 = obj.tiles(obj.pm.adj(ix,2)).tform.T([6]);
%     
%     m1(:,1) = m1(:,1) + x1;
%     m1(:,2) = m1(:,2) + y1;
%     m2(:,1) = m2(:,1) + x2;
%     m2(:,2) = m2(:,2) + y2;
    
    
    for pix = 1:size(m1,1)
        plot(m1, m2, '-');hold on;
    end
end
end