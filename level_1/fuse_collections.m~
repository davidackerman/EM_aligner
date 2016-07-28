function resp = fuse_collections(rcsource, rcfixed, rcmoving, overlap, rcout)
%%% stack fusion
%%% input: fixed collection
%%%        moving collection
%%%        section overlap range
%%%        target output collection (if different from fixed collection)
% clc;clear all;
% 
% % configure source collection
% rcsource.stack          = 'v12_acquire_merged';
% rcsource.owner          ='flyTEM';
% rcsource.project        = 'FAFB00';
% rcsource.service_host   = '10.37.5.60:8080';
% rcsource.baseURL        = ['http://' rcsource.service_host '/render-ws/v1'];
% rcsource.verbose        = 1;
% 
% % configure rough collection
% rcfixed.stack          = ['EXP_dmesh_rough_P1_1_35'];
% rcfixed.owner          ='flyTEM';
% rcfixed.project        = 'test';
% rcfixed.service_host   = '10.37.5.60:8080';
% rcfixed.baseURL        = ['http://' rcfixed.service_host '/render-ws/v1'];
% rcfixed.verbose        = 1;
% rcfixed.nfirst         = 1;
% rcfixed.nlast          = 35;
% 
% % configure rough collection
% rcmoving.stack          = ['EXP_dmesh_rough_P1_20_45'];
% rcmoving.owner          ='flyTEM';
% rcmoving.project        = 'test';
% rcmoving.service_host   = '10.37.5.60:8080';
% rcmoving.baseURL        = ['http://' rcmoving.service_host '/render-ws/v1'];
% rcmoving.verbose        = 1;
% rcmoving.nfirst         = 20;
% rcmoving.nlast          = 45;
% 
% % configure output collection
% rcout.stack          = ['EXP_dmesh_rough_P1_1_45_fused'];
% rcout.owner          ='flyTEM';
% rcout.project        = 'test';
% rcout.service_host   = '10.37.5.60:8080';
% rcout.baseURL        = ['http://' rcmoving.service_host '/render-ws/v1'];
% rcout.verbose        = 1;
% 
% overlap = [20 20];

%% read range
disp('Reading section id ranges...');
[zu1, sID1, sectionId1, z1, ns1] = get_section_ids(rcfixed, overlap(1), overlap(2));
[zu2, sID2, sectionId2, z2, ns2] = get_section_ids(rcmoving, overlap(1), overlap(2));
disp('Done!');
%% find corresponding tile centers
disp('Building list of corresponding tile centers...');
X1 = [];
Y1 = [];
Z1 = [];
X2 = [];
Y2 = [];
Z2 = [];
for zix = 1:numel(zu1)
    %     disp(zix);
    %[x1, y1, tids1] = section_get_tile_centers(rcfixed,zu1(zix));
    %[x2, y2, tids2] = section_get_tile_centers(rcmoving,zu2(zix));
    
    L1 = Msection(rcfixed,zu1(zix));
    L2 = Msection(rcmoving,zu2(zix));
    tids1 = {};
    tids2 = {};
    for tix = 1:numel(L1.tiles), tids1{tix} = L1.tiles(tix).renderer_id;end
    for tix = 1:numel(L2.tiles), tids2{tix} = L2.tiles(tix).renderer_id;end
    [~,ia,ib] = intersect(tids1,tids2);
    
    x1 = L1.X;y1 = L1.Y;
    x2 = L2.X;y2 = L2.Y;
    
    X1 = [X1;x1(ia)];
    Y1 = [Y1;y1(ia)];
    Z1 = [Z1; ones(length(ia),1)];
    X2 = [X2;x2(ib)];
    Y2 = [Y2;y2(ib)];
    Z2 = [Z2; ones(length(ib),1)];
end
disp('Done!');
%% % sosi
% L = Msection(rcfixed, 1);
% rc = rcfixed;
% rc.stack = 'EXP_dmesh_section_one';
% ingest_section_into_renderer_database_overwrite(L,rc, rcsource, pwd, 1);
% Lr= Msection(rc,1);
% delete_renderer_stack(rc);



%% sosi --- visualize
% figure;plot3(X1, Y1, Z1, 'b*');drawnow;pause(1);hold on;plot3(X2,Y2,Z2,'r*');
%% find rotation/translation to overlap tile centers
disp('Determining rotation and translation based on overlap region.....');
A = [X1 Y1];
B = [X2 Y2];
assert(size(A,1)==size(B,1));
centroid_A = mean(A);
centroid_B = mean(B);
N = size(A,1);
H = (A-repmat(centroid_A,N,1))' * (B-repmat(centroid_B,N,1));
[U,S,V] = svd(H);
R = V*U';
if det(R)<0  % detect reflection special case
    V(:,1) = V(:,1) * (-1);
    R = V*U';
end

% find translation
t = -R*centroid_A' + centroid_B';
T = R;
T1 = zeros(3,3);
T1(3,3) = 1;
T1([3 6]) = t(:)';

U = [X1 Y1];% + repmat(t(:)', length(X1),1);
V = [X2 Y2] *R;

x = ones(length(X1),2)\(V-U);
V = V - repmat(x(1,:), length(X1),1);

T2 = zeros(3,3);
T2(3,3) = 1;
T2([3 6]) = x(1,:);

disp('Done!');

%% %% do again
% B = [U(:,1) U(:,2)];
% A = [V(:,1) V(:,2)];
% assert(size(A,1)==size(B,1));
% centroid_A = mean(A);
% centroid_B = mean(B);
% N = size(A,1);
% H = (A-repmat(centroid_A,N,1))' * (B-repmat(centroid_B,N,1));
% [U,S,V] = svd(H);
% R2 = V*U';
% if det(R)<0  % detect reflection special case
%     V(:,1) = V(:,1) * (-1);
%     R2 = V*U';
% end
% 
% % find translation
% t = -R2*centroid_A' + centroid_B';
% T = R2;
% T(3,3) = [1];
% T([3 6]) = [t(1) t(2)];
% 
% U = A *R2 + repmat(t(:)', length(X1),1);
% V = B;
% % 
% x = ones(length(X1),2)\(U-V);
% U = U - repmat(x(1,:), length(X1),1);


% figure;plot3(U(:,1), U(:,2), Z1, 'b*');drawnow;pause(1);hold on;plot3(V(:,1), V(:,2) ,Z2,'r*');drawnow;
% 
% figure;L12 = Msection(rcfixed, 20);
% show_map(L12);
% hold on;plot3(V(:,1), V(:,2) ,Z2,'r*');drawnow;
%% read collections --- L11 L12 L21 L22


% [zu1, sID1, sectionId1, z1, ns1] = get_section_ids(rcfixed, rcfixed.nfirst, overlap(1));
% tiles11 = [];
% for ix = 1:numel(zu1)-1
%     L11 = Msection(rcfixed, zu1(ix));
%     tiles11 = [tiles11 L11.tiles];
% end
% L11 = Msection(tiles11);    % part that does not get modified


disp('Reading tile collections ...');
disp('Overlap region of fixed collection.');
[zu1, sID1, sectionId1, z1, ns1] = get_section_ids(rcfixed, overlap(1), overlap(2));
tiles12 = [];
for ix = 1:numel(zu1)
    L12 = Msection(rcfixed, zu1(ix));
    tiles12 = [tiles12 L12.tiles];
end
L12 = Msection(tiles12);  % part that overlaps from fixed

disp('Overlap region of moving collection.');
[zu2, sID2, sectionId2, z2, ns2] = get_section_ids(rcmoving, rcmoving.nfirst, overlap(2));
tiles21 = [];
for ix = 1:numel(zu2)
    L21 = Msection(rcmoving, zu2(ix));    
    tiles21 = [tiles21 L21.tiles];
end
L21 = Msection(tiles21); % part that overlaps from moving (gets rotated and translated)

disp('Non-overlap region of moving collection.');
[zu2, sID2, sectionId2, z2, ns2] = get_section_ids(rcmoving, overlap(2)+1, rcmoving.nlast);
tiles22 = [];
for ix = 1:numel(zu2)
    L22 = Msection(rcmoving, zu2(ix));    
    tiles22 = [tiles22 L22.tiles];
end
L22 = Msection(tiles22); % part that does not ovelap from moving (but gets rotated and translated)
disp('Done!');
%% transform all rcmoving
disp('Applying transformations to moving tiles...');
tiles21t = tiles21;
tiles22t = tiles22;
H = tiles21(1).H;
W = tiles21(1).W;

Tr = zeros(3,3);Tr(3,3) = 1;
Tr(1:2,1:2) = R;
delta = [W H];
%delta = [0 0];
T3 = T2;
T3([3 6]) = T3([3 6]) + delta;



T = -T3  + Tr;
T(3,3) = 1;

for tix = 1:numel(L21.tiles)
    tiles21t(tix).tform.T = tiles21(tix).tform.T * T;
end

for tix = 1:numel(L22.tiles)
    tiles22t(tix).tform.T = tiles22(tix).tform.T * T;
end
L21t = Msection(tiles21t);
L22t = Msection(tiles22t);

disp('Done!');
%% interpolate: tiles in L21t need to be modified to comply with interpolation
% to generate a new set Linterp to be used
% tids1 = {};
% tids2 = {};
% for tix = 1:numel(tiles12), tids1{tix} = tiles12(tix).renderer_id;end
% for tix = 1:numel(tiles21t), tids2{tix} = tiles21t(tix).renderer_id;end
% [~,ia,ib] = intersect(tids1,tids2);
% tiles12 = tiles12(ia);
% tiles21t = tiles21t(ib);



disp('Interpolating tile transformations within overlap region...');
% % determine intersecting tiles for each section and interpolate
% [zu1, sID1, sectionId1, z1, ns1] = get_section_ids(rcfixed, overlap(1), overlap(2));
[zu2, sID2, sectionId2, z2, ns2] = get_section_ids(rcmoving, rcmoving.nfirst, overlap(2));
dlambda = 1/numel(zu2);
lambda = 0;
del_ix = [];
for zix = 1:numel(zu2)
    disp([zix zu2(zix) lambda]);
    zcurr = zu2(zix);
    tindx = find([tiles21t(:).z]==zcurr);
    for tix = 1:numel(tindx)
        if isKey(L12.map_renderer_id, tiles21t(tindx(tix)).renderer_id)
            ind1 = L12.map_renderer_id(tiles21t(tindx(tix)).renderer_id);
            if zcurr==26, disp([tiles21t(tindx(tix)).tform.T(:)' - tiles12(ind1).tform.T(:)']);end
            tiles21t(tindx(tix)).tform.T(1:3,1:2) = tiles12(ind1).tform.T(1:3,1:2).* (1-lambda) +...
                tiles21t(tindx(tix)).tform.T(1:3,1:2).* lambda;
        else
            del_ix = [del_ix tindx(tix)];
        end
    end
    lambda = lambda + dlambda;
end

L21t = Msection(tiles21t);
L21t.tiles(del_ix) = [];
disp('Done!');
%%
% L = Msection([L11.tiles L21t.tiles L22t.tiles]);
disp('Assembling slab composed of transformed (and interpolated) sections ...');
L = Msection([L21t.tiles L22t.tiles]);
disp('Done!');
%%
% close all;
% %%% result compare overlap-1 and overlap
% zmL = split_z(L);
% 
% show_map(zmL(20));


% show_map(zmL(overlap(1)-1));
% hold on; drawnow
% show_map(zmL(overlap(1)));
% 
% %%% compare original overlap-1 with original overlap in rcfixed
% figure;
% zmL11 = split_z(L11);
% show_map(zmL11(overlap(1)-1));
% hold on;
% L12 = Msection(rcfixed, overlap(1));
% show_map(L12);hold on;plot3(V(:,1), V(:,2) ,Z2,'r*');drawnow;
% 
% %%% compare result overlap(1) with original overlap(1) (in rcfixed)
% figure;
% zmL21t = split_z(L21t);
% show_map(zmL21t(1));
% hold on;
% L12 = Msection(rcfixed, overlap(1));
% show_map(L12);drawnow;



%% export to collection
%%% ingest into renderer database as rough collection
disp('Ingesting into renderer database using append.');
resp = ingest_section_into_LOADING_collection(L,rcout, rcsource, pwd);
resp = set_renderer_stack_state_complete(rcout);
disp('Done!');
% 

% resp = set_renderer_stack_state_complete(rcout);
% %%
% L21t = Msection(tiles21t);
% L22t = Msection(tiles22t);
% tids1 = {};
% tids2 = {};
% for tix = 1:numel(L12.tiles), tids1{tix} = L12.tiles(tix).renderer_id;end
% for tix = 1:numel(L21.tiles), tids2{tix} = L21.tiles(tix).renderer_id;end
% [~,ia,ib] = intersect(tids1,tids2);
% x1 = L12.X;
% y1 = L12.Y;
% x2 = L21.X;
% y2 = L21.Y;
% 
% X1 = x1(ia);
% Y1 = y1(ia);
% X2 = x2(ib);
% Y2 = y2(ib);
% 
% 
% plot(X1, Y1, 'b*');drawnow;pause(1);hold on;plot(X2,Y2,'r*');


















