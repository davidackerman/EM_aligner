function resp = fuse_collections(rcsource, rcfixed, rcmoving, overlap, rcout, collection_start)
%%% stack fusion
%%% input: fixed collection
%%%        moving collection
%%%        section overlap range
%%%        target output collection (if different from fixed collection)


%% read range
disp('Reading section id ranges...');
[zu1, sID1, sectionId1, z1, ns1] = get_section_ids(rcfixed, overlap(1), overlap(2));
[zu2, sID2, sectionId2, z2, ns2] = get_section_ids(rcmoving, overlap(1), overlap(2));
disp('Done!');
%% find corresponding tile centers
disp('Building list of corresponding tile centers...');
tic
X1 = {};
Y1 = {};
X2 = {};
Y2 = {};

parfor zix = numel(zu1)    
    %disp([zu1(zix) zu2(zix)]);
    %[x1, y1, tids1] = section_get_tile_centers(rcfixed,zu1(zix));
    %[x2, y2, tids2] = section_get_tile_centers(rcmoving,zu2(zix));

%     L1 = Msection(rcfixed,zu1(zix));
%     tids1 = {};
%     for jix = 1:numel(L1.tiles)
%         tids1{jix} = L1.tiles(jix).renderer_id;
%     end
%     x1 = L1.X;y1 = L1.Y;
%     L2 = Msection(rcmoving,zu2(zix));
%     tids2 = {};
%     for jix = 1:numel(L2.tiles)
%         tids2{jix} = L2.tiles(jix).renderer_id;
%     end
%     x2 = L2.X;y2 = L2.Y;
% 
%     [x1, y1, tids1, L1f] = get_tile_centers(rcfixed, zu1(zix));
%     disp([L1.X x1 L1.Y y1]);
    
    [x1, y1, tids1] = get_tile_centers(rcfixed, zu1(zix));
    [x2, y2, tids2] = get_tile_centers(rcmoving, zu2(zix));

    [~,ia,ib] = intersect(tids1,tids2);
    
        X1{zix} = [x1(ia)'];
        Y1{zix} = [y1(ia)'];
        X2{zix} = [x2(ib)'];
        Y2{zix} = [y2(ib)'];
end
X1 = cell2mat(X1)';
Y1 = cell2mat(Y1)';
X2 = cell2mat(X2)';
Y2 = cell2mat(Y2)';

toc
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
disp('-- Building linear system for rotation');
N = size(A,1);
H = (A-repmat(centroid_A,N,1))' * (B-repmat(centroid_B,N,1));
[U,S,V] = svd(H);
R = V*U';
if det(R)<0  % detect reflection special case
    V(:,1) = V(:,1) * (-1);
    R = V*U';
end
disp('-- building linear system for translation');
% find translation
% t = -R*centroid_A' + centroid_B';
% T = R;
% T2 = zeros(3,3);
% T2(3,3) = 1;
% T2([3 6]) = -t(:)';

U = [X1 Y1];% + repmat(t(:)', length(X1),1);
V = [X2 Y2] *R;
warning off; x = ones(length(X1),2)\(V-U);warning on;
V = V - repmat(x(1,:), length(X1),1);
T2 = zeros(3,3);
T2(3,3) = 1;
T2([3 6]) = x(1,:);

disp('Done!');
%% inspect result:

% % %
% figure;
A = [X1 Y1];
B = [X2 Y2];

plot(A(:,1), A(:,2), '*b');
hold on;
[tform] = cp2tform(B, A, 'nonreflective similarity');
[Btx, Bty] = tformfwd(tform, B(:,1), B(:,2));
plot(Btx, Bty, '*r');
%%

% figure;
 A = [X1 Y1];
 B = [X2 Y2];
% 
% plot(A(:,1), A(:,2), '*b');
% hold on;
 [tform] = cp2tform(B, A, 'nonreflective similarity');
% Btr = [B ones(size(B,1),1)]*tform.tdata.T(1:3,1:3) ;% + repmat(tform.tdata.T([3 6]), size(B,1),1);
% plot(Btr(:,1), Btr(:,2), '*r');
% [Btx Btr(:,1)]

T = tform.tdata.T;

%% read collections --- (L11) L12 L21 L22

tic
disp('Reading tile collections ...');

if collection_start
disp('-- Non-overlap region of fixed collection.');
[~, tiles11] = get_slab_tiles(rcfixed, rcfixed.nfirst, overlap(1)-1);
end


disp('-- Overlap region of fixed collection.');
[~, tiles12] = get_slab_tiles(rcfixed, overlap(1), overlap(2));


disp('-- Overlap region of moving collection.');
[~, tiles21] = get_slab_tiles(rcmoving, overlap(1), overlap(2));

disp('-- Non-overlap region of moving collection.');
[~, tiles22] = get_slab_tiles(rcmoving, overlap(2)+1, rcmoving.nlast);

disp('Done!');
toc

%% transform all rcmoving
tic
disp('Applying transformations to moving tiles...');
% initialize transformed tile arrays
tiles21t = tiles21;
tiles22t = tiles22;

% H = tiles21(1).H;
% W = tiles21(1).W;
% 
% % assemble full trnsformation matrix rotation and translation
% Tr = zeros(3,3);Tr(3,3) = 1;
% Tr(1:2,1:2) = R;
% delta = [W H];
% T3 = T2;
% T3([3 6]) = T3([3 6]) + delta;
% T = -T3  + Tr;
% T(3,3) = 1;

% apply transformations
for tix = 1:numel(tiles21)
    tiles21t(tix).tform.T = tiles21(tix).tform.T *T;

end

for tix = 1:numel(tiles22)
    tiles22t(tix).tform.T = tiles22(tix).tform.T *T;

end
L21t = Msection(tiles21t);
L22t = Msection(tiles22t);
disp('Done!');
toc




%% interpolate: tiles in L21t need to be modified to comply with interpolation



tic
disp('Interpolating tile transformations within overlap region...');
% % determine intersecting tiles for each section and interpolate
% [zu1, sID1, sectionId1, z1, ns1] = get_section_ids(rcfixed, overlap(1), overlap(2));
[zu2, sID2, sectionId2, z2, ns2] = get_section_ids(rcmoving, rcmoving.nfirst, overlap(2));
dlambda = 1/numel(zu2);
lambda = 0;
del_ix = [];
L12 = Msection(tiles12);
for zix = 1:numel(zu2)
    %disp([zix zu2(zix) lambda]);
    zcurr = zu2(zix);
    tindx = find([tiles21t(:).z]==zcurr);
    for tix = 1:numel(tindx)
        if isKey(L12.map_renderer_id, tiles21t(tindx(tix)).renderer_id)
            ind1 = L12.map_renderer_id(tiles21t(tindx(tix)).renderer_id);
            %if zcurr==26, disp([tiles21t(tindx(tix)).tform.T(:)' - tiles12(ind1).tform.T(:)']);end
            tiles21t(tindx(tix)).tform.T(1:3,1:2) = tiles12(ind1).tform.T(1:3,1:2).* (1-lambda) +...
                tiles21t(tindx(tix)).tform.T(1:3,1:2).* lambda;
        else
            del_ix = [del_ix tindx(tix)];
        end
    end
    lambda = lambda + dlambda;
end

tiles21t(del_ix) = [];

disp('Done!');
toc
%%
disp('Assembling slab composed of transformed (and interpolated) sections ...');
tic
if collection_start,
         
    disp('-- collection start: assembling full set');
    L = Msection([tiles11(:)' tiles21t(:)' tiles22t(:)']);
    %L = Msection([L11.tiles; L21t.tiles; L22t.tiles]);
    
else
    L = Msection([tiles21t(:)' tiles22t(:)']);
end
toc
disp('Done!');

%%
% close all;
% %%% result compare overlap-1 and overlap
% zmL = split_z(L);
% show_map(zmL(19));
% drawnow;pause(1);hold on;
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
%%% compare result overlap(1) with original overlap(1) (in rcfixed)
figure;
zmL21t = split_z(L21t);
show_map(zmL21t(1));
hold on;
L12 = Msection(rcfixed, overlap(1));
show_map(L12);drawnow;



%% export to collection
%%% ingest into renderer database as rough collection
tic
disp('Ingesting into renderer database using append.');
%resp = ingest_section_into_LOADING_collection(L,rcout, rcsource, pwd);


% disp('-- Translate to positive space');
% L = translate_to_origin(L);
disp('-- Splitting into sections to prepare for distributed ingestion');
zmL = split_z(L);
disp('-- Start distributed process to populate new renderer collection');
resp_append = {};
translate_to_positive_space = 0;
complete = 0;
parfor ix = 1:numel(zmL)
    ingest_section_into_renderer_database(zmL(ix), ...
        rcout, rcsource, pwd, translate_to_positive_space, complete);
end


resp = set_renderer_stack_state_complete(rcout);
disp('Done!');
toc
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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%













