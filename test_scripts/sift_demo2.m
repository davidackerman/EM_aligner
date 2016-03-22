
% t1 = L(l1(pix)).tiles(id1(pix));
% t2 = L(l2(pix)).tiles(id2(pix));
% 
% t1.fetch_local = 0;
% t2.fetch_local = 0;

% im1 = get_image(t1);
% im2 = get_image(t2);


 t1 = L(1).tiles(10);
 t2 = L(1).tiles(11);
% 
% im1 = imread(t1.path);
% im2 = imread(t2.path);

im1 = get_image(t1);
im2 = get_image(t2);

scale = 0.25;
thresh = 5;
indx1 = find(im1<thresh);r = rand(size(indx1));im1(indx1) = r * double(max(im1(:)));
indx2 = find(im2<thresh);r = rand(size(indx2));im2(indx2) = r * double(max(im2(:)));
%figure;imshowpair(im1, im2, 'montage');

I1=imresize(double(im1), scale) ; % I1=I1(1:2:end,:) ;
I2=imresize(double(im2), scale) ; % I2=I2(1:2:end,:) ;
I1c=double(im1)/255.0 ;
I2c=double(im2)/255.0 ;


I1=imsmooth(I1,.1) ;
I2=imsmooth(I2,.1) ;

I1=I1-min(I1(:)) ;
I1=I1/max(I1(:)) ;
I2=I2-min(I2(:)) ;
I2=I2/max(I2(:)) ;

S=8 ;

fprintf('Computing frames and descriptors.\n') ;
[frames1,descr1,gss1,dogss1] = sift( I1, 'Verbosity', 1, 'Threshold', ...
                                     0.005, 'NumLevels', S ) ;
[frames2,descr2,gss2,dogss2] = sift( I2, 'Verbosity', 1, 'Threshold', ...
                                     0.005, 'NumLevels', S ) ;

% figure(11) ; clf ; plotss(dogss1) ; colormap gray ;
% figure(12) ; clf ; plotss(dogss2) ; colormap gray ;
% drawnow ;
% 
% figure(2) ; clf ;
% tightsubplot(1,2,1) ; imagesc(I1) ; colormap gray ; axis image ;
% hold on ;
% h=plotsiftframe( frames1 ) ; set(h,'LineWidth',2,'Color','g') ;
% h=plotsiftframe( frames1 ) ; set(h,'LineWidth',1,'Color','k') ;
% 
% tightsubplot(1,2,2) ; imagesc(I2) ; colormap gray ; axis image ;
% hold on ;
% h=plotsiftframe( frames2 ) ; set(h,'LineWidth',2,'Color','g') ;
% h=plotsiftframe( frames2 ) ; set(h,'LineWidth',1,'Color','k') ;

fprintf('Computing matches.\n') ;
% By passing to integers we greatly enhance the matching speed (we use
% the scale factor 512 as Lowe's, but it could be greater without
% overflow)
descr1=uint8(512*descr1) ;
descr2=uint8(512*descr2) ;
tic ;
matches=siftmatch( descr1, descr2, 2.0) ;
fprintf('Matched in %.3f s\n', toc) ;

figure(3) ; clf ;
plotmatches(I1c,I2c,frames1(1:2,:),frames2(1:2,:),matches,...
  'Stacking','v') ;
drawnow ;

% % Movie
% figure(4) ; set(gcf,'Position',[10 10 1024 512]) ;
% figure(4) ; clf ;
% tightsubplot(1,1);
% imagesc(I1) ; colormap gray ; axis image ; hold on ;
% h=plotsiftframe( frames1 ) ; set(h,'LineWidth',1,'Color','g') ;
% h=plot(frames1(1,:),frames1(2,:),'r.') ;
% MOV(1)=getframe ;
% 
% figure(4) ; clf ;
% tightsubplot(1,1);
% imagesc(I2) ; colormap gray ; axis image ; hold on ;
% h=plotsiftframe( frames2 ) ; set(h,'LineWidth',1,'Color','g') ;
% h=plot(frames2(1,:),frames2(2,:),'r.') ;
% MOV(2)=getframe ;
