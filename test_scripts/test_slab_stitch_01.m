%% The full stitching pipeline working in conjunction with the Renderer database and the point-match database
%% for stitching a three-layer slab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [0] configure collections and prepare quantities

% configure source collection
rcsource.stack = 'v12_acquire_merged';
rcsource.owner='flyTEM';
rcsource.project= 'FAFB00';
rcsource.server = 'http://10.37.5.60:8080/render-ws/v1';
rcsource.service_host = '10.37.5.60:8080';

rcsource.verbose = 1;
rcsource.baseURL        = rcsource.server;
rcsource.source_stack   = rcsource.stack;
rcsource.verbose        = 1;

% configure montage target collection
rctarget_montage = rcsource;
rctarget_montage.project = 'test';

rctarget_rough = rcsource;
rctarget_rough.project = 'test';


% configure point-match collection
pm.server = 'http://10.40.3.162:8080/render-ws/v1';
pm.owner  = 'flyTEM';
pm.match_collection = 'v12_SURF';

rough_align_method = 'SIFT'; % SIFT uses Stephan's code via Trautman's scirpt. The other option is to use SURF (matlab code)
rough_align_scale  = 0.05;
generate_rough_alignment_preview = 1;

nfirst = 1245;
nlast  = 1247;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% get the list of zvalues and section ids within the z range between nfirst and nlast (inclusive)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
urlChar = sprintf('%s/owner/%s/project/%s/stack/%s/sectionData', ...
     rcsource.server, rcsource.owner, rcsource.project, rcsource.stack);
 j = webread(urlChar);
 sectionId = {j(:).sectionId};
 z         = [j(:).z];
 indx = find(z>=nfirst & z<=nlast);
 
 sectionId = sectionId(indx);% determine the sectionId list we will work with
 z         = z(indx);        % determine the zvalues (this is also the spatial order)
%% [1] generate montage for individual sections and generate montage collection
L = Msection;
L(numel(z)) = Msection;
for lix = 1:numel(z)
    L(lix) = Msection(rcsource, z(lix));
    L(lix).dthresh_factor = 3;
    L(lix) = update_XY(L(lix));
    L(lix)  = update_adjacency(L(lix));
    [L(lix), js] = alignTEM_inlayer(L(lix));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% ingest js into point matches database
    fn = ['temp_' num2str(randi(100000)) '_' num2str(lix) '.json'];
    fid = fopen(fn, 'w');
    fwrite(fid,js);
    fclose(fid);
    urlChar = sprintf('%s/owner/%s/matchCollection/%s/matches/', ...
        pm.server, pm.owner, pm.match_collection);
    cmd = sprintf('curl -X PUT --connect-timeout 30 --header "Content-Type: application/json" --header "Accept: application/json" -d "@%s" "%s"',...
        fn, urlChar);
    [a, resp]= evalc('system(cmd)');
    delete(fn);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
opts.outlier_lambda = 1e3;  % large numbers result in fewer tiles excluded
mL = concatenate_tiles(L, opts.outlier_lambda);
collection = ['EXP_v12_SURF_montage_' num2str(nfirst) '_' num2str(nlast)];
ingest_section_into_renderer_database(mL, rctarget_montage, rcsource, pwd);
mL = update_tile_sources(mL, rctarget_montage.owner, rctarget_montage.project, collection, rcsource.server);
L = split_z(mL);
%% [2] generate montage-scapes and montage-scape point-matches 
rc_montage = rcsource;
rc_montage.stack = collection;
if strcmp('SIFT', rough_align_method)
    

    script_full_path = '/tier2/flyTEM/khairy/sw_work/spark_montage_scape_pm_generation/run_FAFB_montage_scapes_experiments.sh';
    sp.fd_size = 4;
    sp.min_sift_scale = 0.75;
    sp.max_sift_scale = 1.0;
    sp.steps = 3;
    scale = 0.1;
    similarity_range = 6;
    dir_output = '/nobackup/flyTEM/spark_montage';    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% run script by Eric T.
    %%% separately for now
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp( '-------------------- run manually to completion ... then press any key');pause(0);
    fn_matches = '/nobackup/flyTEM/spark_montage/test/EXP_v12_SURF_montage_1245_1247/scale_0.05_sift_8_0.55_1.0_3/solver_1245.00_1247.00/matches.txt';
    if exist(fn_matches, 'file'), delete(fn_matches);end % remove the file it it already exists
    cmd_str = ...
        sprintf('/tier2/flyTEM/khairy/sw_work/spark_montage_scape_pm_generation/run_FAFB_montage_scapes_experiments_v12_1245_1247.sh 10.40.3.162:8080');
    [a, resp_str] = system(cmd_str);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    while exist(fn_matches,'file')~=2
        disp('Waiting for files to finish generating');
        pause(30);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    L2 = get_montage_scapes_SIFT_point_matches(...
        script_full_path, rc_montage,...
        [1 nlast-nfirst+1], sp, scale, similarity_range, dir_output);
    for tix = 1:numel(L2.tiles)
        L2.tiles(tix).fetch_local = 1;
    end
%     fac = 4;
%     for pmix = 1:size(L2.pm.adj,1)
%         L2.pm.M{pmix,1} = L2.pm.M{pmix,1} * fac;
%         L2.pm.M{pmix,2} = L2.pm.M{pmix,2} * fac;
%     end
    %%% filter point matches using RANSAC
    geoTransformEst = vision.GeometricTransformEstimator; % defaults to RANSAC
    geoTransformEst.Method = 'Random Sample Consensus (RANSAC)';%'Least Median of Squares';
    geoTransformEst.Transform = 'Affine';%'Nonreflective similarity';%'Affine';%
    geoTransformEst.NumRandomSamplingsMethod = 'Desired confidence';
    geoTransformEst.MaximumRandomSamples = 1000;
    geoTransformEst.DesiredConfidence = 99.95;
    for pmix = 1:size(L2.pm.M,1)
        m1 = L2.pm.M{pmix,1};
        m2 = L2.pm.M{pmix,2};
        % Invoke the step() method on the geoTransformEst object to compute the
        % transformation from the |distorted| to the |original| image. You
        % may see varying results of the transformation matrix computation because
        % of the random sampling employed by the RANSAC algorithm.
        [tform_matrix, inlierIdx] = step(geoTransformEst, m2, m1);
        m1 = m1(inlierIdx,:);
        m2 = m2(inlierIdx,:);
        L2.pm.M{pmix,1} = m1;
        L2.pm.M{pmix,2} = m2;
        w = L2.pm.W{pmix};
        L2.pm.W{pmix} = w(inlierIdx);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
else
    disp('Using Matlab''s SURF to determine rough alignment');
    maxTiles = 50;
    im = {};
    fn = {};
    dir_tmp = '/nobackup/flyTEM/khairy/temp_stitch';
    kk_mkdir(dir_tmp);
    ms_tiles = tile;
    ms_tiles(numel(z)) = tile;
    for lix = 1:numel(z)
        %%% get montage scape from renderer service
        urlChar = sprintf('%s/owner/%s/project/%s/stack/%s/z/%s/jpeg-image?scale=%s&filter=%s&maxTileSpecsToRender=%s', ...
            rc_montage.server, rc_montage.owner, rc_montage.project, rc_montage.stack, num2str(z(lix)), num2str(rough_align_scale), 'true', num2str(maxTiles));
        im1 = webread(urlChar);
        im1 = rgb2gray(im1);
        indx = find(im1==0);
        im1(indx) = uint8(rand(numel(indx),1))*max(im1(:));
        im{lix} = im1;
        %%% persist to disc
        fn{lix} = [dir_tmp '/montage_scape_' num2str(lix) '.jpg'];
        imwrite(im{lix}, fn{lix});
    end
    % % Generate montage features and collect in an Msection object
    for lix = 1:numel(z)
        %%% make tile out of montage scape image
        t = tile;
        t.z = lix;
        t.id = lix * 1000000000 + randi(100000);
        t.path = fn{lix};
        t.fetch_local = 1;
        t.rot = 0;
        t.featuresMethod = 'SURF';
        t.SURF_NumOctaves = 2;
        t.SURF_NumScaleLevels = 3;
        t.SURF_MetricThreshold = 500;
        
        t = get_features(t, 3);
        ms_tiles(lix) = t;
    end
    msL = Msection(ms_tiles);
    % generate point matches
    min_pm = 1;
    [L2] = generate_point_matches(msL, min_pm);
end


% check the point matches
for lix = 1%:size(L2.pm.adj,1)
    ix1 = L2.pm.adj(lix,1);
    ix2 = L2.pm.adj(lix,2);
    t1 = L2.tiles(ix1);
    t2 = L2.tiles(ix2);
    M = L2.pm.M(lix,:);
    show_feature_point_correspondence(t1,t2,M);title([num2str(ix1) '   ' num2str(ix2)]);
    drawnow;
end
%% [3] rough alignment of montage-scapes

% solve
[mLR, errR, mLS] = get_rigid_approximation(L2);
[mL3, errA] = solve_affine_explicit_region(mLR);
if generate_rough_alignment_preview
    mL3.update_tile_info_switch = 1;
    mL3 = update_tile_info(mL3);
% render to look at quality of alignment of montage scapes
mL3 = update_XY(mL3);
mL3 = get_bounding_box(mL3);
Wbox = [mL3.box(1) mL3.box(3) mL3.box(2)-mL3.box(1) mL3.box(4)-mL3.box(3)];
disp(Wbox);
dir_out = ['/nobackup/flyTEM/khairy/temp_stitch/render_test_v9_mL3_' num2str(nfirst) '_to_' num2str(nlast)];
kk_mkdir(dir_out);
disp(dir_out);
kk_clock();
disp('Rendering ...');tic;
[im, Wbox, imlabel, M] = render_poly_07(mL3.tiles, 1.0, Wbox, 0,[], dir_out);
disp('Done!');toc;
end
mL3s = split_z(mL3);
%% [4] apply rough alignment to sections and generate "rough_aligned" collection  %% %%%%%% sosi
for lix = 1:numel(L), L(lix) = get_bounding_box(L(lix));end
Wbox = [mL3.box(1) mL3.box(3) mL3.box(2)-mL3.box(1) mL3.box(4)-mL3.box(3)];disp(Wbox);
wb1 = Wbox(1);
wb2 = Wbox(2);
L3 = L;
fac = 0.25;

smx = [fac 0 0; 0 fac 0; 0 0 1]; %scale matrix
invsmx = [1/fac 0 0; 0 1/fac 0; 0 0 1];
tmx2 = [1 0 0; 0 1 0; -wb1 -wb2 1]; % translation matrix for montage_scape stack

for lix = 1:numel(L)
    b1 = L(1).box;
    dx = b1(1);dy = b1(3);
    tmx1 = [1 0 0; 0 1 0; -dx -dy 1];  % translation matrix for section box
    for tix = 1:numel(L3(lix).tiles)
        newT = L3(lix).tiles(tix).tform.T * tmx1 * smx * mL3s(lix).tiles(1).tform.T * tmx2 * (invsmx);
        L3(lix).tiles(tix).tform.T = newT;     
        
% %         %%%%%% also transform locations of features -------- sosi
% %         featurePointT = tmx1 * smx * mL3s(lix).tiles(1).tform.T * tmx2;
% %         validPoints.Location = ...
% %             [L3(lix).tiles(tix).validPoints.Location ones(length(L3(lix).tiles(tix).validPoints.Location),1)]...
% %             * featurePointT;
% %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    L3(lix) = update_XY(L3(lix));
end

%%%%%%%%%%% sosi
% figure;show_map(L3(4));hold on; drawnow;pause(2); show_map(L3(5));drawnow;
%figure;show_map(L3(3));hold on; drawnow;pause(2); show_map(L3(5));drawnow;
% L3(4).tiles(9).tform.T
% L3(5).tiles(7).tform.T
% L(5).tiles(7).tform.T
% 
% mL3s(4).tiles(1).tform.T
%mL3s(5).tiles(1).tform.T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opts.outlier_lambda = 1e3;  % large numbers result in fewer tiles excluded
mL = concatenate_tiles(L3, opts.outlier_lambda);
collection = ['EXP_v12_SURF_rough_align' num2str(nfirst) '_' num2str(nlast)];
ingest_section_into_renderer_database(mL,rctarget_rough, rcsource, pwd);
mL = update_tile_sources(mL, rctarget_rough.owner, rctarget_rough.project, collection, rcsource.server);
L3 = split_z(mL);
%% [5] Determine potential tile partners crosslayer (determine the list of cross-layer tile partners)
% first determine the list of section pairs
dthresh_factor = 4;
depth = 1;
top = numel(L3);
bottom = 1;
cs = [];  % compare section list
counter = 1;
for uix = top:-1:bottom
    lowest = uix-depth;
    if lowest<bottom, lowest = bottom;end
    for vix = uix:-1:lowest
        if vix<uix
            cs(counter,:) = [uix vix];
            counter = counter + 1;
        end
    end
end
disp(cs);

% loop over cs pairs and determine potential tile pairs
l1 = [];
l2 = [];
id1 = [];
id2 = [];
rid1 = {};
rid2 = {};
H = L3(1).tiles(1).H;
W = L3(1).tiles(1).W;
for dix = 1:size(cs,1)
    lix1 = cs(dix,1);
    lix2 = cs(dix,2);
    a = [L3(lix1).X L3(lix1).Y];
    b = [L3(lix2).X L3(lix2).Y];
    d = pdist2(a,b);        % depends on statistics toolbox  -------- Sosi: not good for large numbers of tiles
    dthresh = sqrt(H^2 + W^2) * dthresh_factor;   % diagonal of the tile times factor
    A = sparse(triu(d<dthresh,0)); % generate adjacency matrix
    [r c] = ind2sub(size(A), find(A));
    l1 = [l1 lix1 * ones(1,size(r(:),1))];
    l2 = [l2 lix2 * ones(1,size(c(:),1))];
    id1 = [id1 r(:)'];
    id2 = [id2 c(:)'];
    rid1 = [rid1 {L3(lix1).tiles(r).renderer_id}];
    rid2 = [rid2 {L3(lix2).tiles(c).renderer_id}];
end
%% [6] fine alignment (cross-layer point-matches)
% for lix = 1:numel(L3)
%     for tix = 1:numel(L3(lix).tiles)
%         L3(lix).tiles(tix).fetch_local = 0;
%         L3(lix).tiles(tix).SURF_NumOctaves = 2;
%         L3(lix).tiles(tix).SURF_NumScaleLevels = 8;
%         L3(lix).tiles(tix).SURF_MetricThreshold= 500;
%     end
%     L3(lix) = calculate_tile_features(L3(lix), 'true', 1, 5000);
% end

minpm = 1;
M = cell(numel(id1),2);
adj = zeros(numel(id1),2);
W = cell(numel(id1),1);
np = zeros(numel(id1),1);
delpix = zeros(numel(id1),1, 'uint32');
parfor_progress(numel(id1));
parfor pix = 1: numel(id1)
    disp([pix l1(pix) l2(pix) id1(pix) id2(pix)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t1 = L3(l1(pix)).tiles(id1(pix));
    t2 = L3(l2(pix)).tiles(id2(pix)); 
    t1.fetch_local = 0;
    t2.fetch_local = 0;
   
    im1 = get_image(t1);
    im2 = get_image(t2);
    [m1, m2, w, imt1] = kk_dftregistration(im1, im2);
    
    %%%%%% sosi
%      figure; showMatchedFeatures(im1, im2, m1, m2, 'montage');
%      figure;imshowpair(im1, imt1, 'montage');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    M(pix,:) = {[m1],[m2]};
    
    
    %%%%%%%% sosi
%     show_feature_point_correspondence(t1,t2,M(pix,:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    adj(pix,:) = [id1(pix) id2(pix)];
    W(pix) = {[ones(size(m1,1),1) * 1/(1+ abs(L3(l1(pix)).tiles(id1(pix)).z-L3(l2(pix)).tiles(id2(pix)).z))]};
    np(pix)  = size(m1,1);
    %%%% mark for removal point-matches that don't have enough point pairs
    if size(m1,1)<minpm
        delpix(pix) = 1;
    end
    parfor_progress;
end
parfor_progress(0);
if isempty(M), error('No matches found');end;
disp('Done!');


delpix = logical(delpix);
M(delpix,:) = [];
adj(delpix,:) = [];
W(delpix) = [];
np(delpix) = [];

l1(delpix) = [];
l2(delpix) = [];
id1(delpix) = [];
id2(delpix) = [];
rid1(delpix) = [];
rid2(delpix) = [];


pm.M = M;
pm.adj = adj;
pm.W = W;
pm.np = np;

%%%%%%%%%% sosi
% % %% check the point matches
% for lix = 1:size(pm.adj,1)
%     ix1 = pm.adj(lix,1);
%     ix2 = pm.adj(lix,2);
%     t1 = L3(l1(lix)).tiles(ix1);
%     t2 = L3(l2(lix)).tiles(ix2);
%     M = pm.M(lix,:);
%     str = [num2str(size(pm.M{lix,1},1)), 'Section: ' t1.sectionId ' tile: ' num2str(ix1) '   ' 'Section: ' t2.sectionId ' tile: ' num2str(ix2)];
% %     close all;    show_feature_point_correspondence(t1,t2,M);
% %     title(str);    %truesize(gcf, [1000 1000]);
% %     drawnow;pause(3);
% disp(str);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [7] transform point matches to lens-corrected-only stage and ingest into pm database
for pix = 1:size(pm.adj,1)
    t1 = L3(l1(pix)).tiles(id1(pix));
    t2 = L3(l2(pix)).tiles(id2(pix)); 
    pm1 = [pm.M{pix,1} ones(size(pm.M{pix,1},1),1)] * inv(t1.tform.T);
    pm2 = [pm.M{pix,2} ones(size(pm.M{pix,2},1),1)] * inv(t2.tform.T);
    
    pm.M{pix,1} = pm1(:,1:2);
    pm.M{pix,2} = pm2(:,1:2);
end
% generate json point-match data
counter = 1;
M = pm.M;
adj = pm.adj;
for mix = 1:size(M,1)
    indx1 = adj(mix,1);
    indx2 = adj(mix,2);
    tid1 = [L3(l1(mix)).tiles(indx1).renderer_id];
    tid2 = [L3(l2(mix)).tiles(indx2).renderer_id];

    MP{counter}.pz = L3(l1(mix)).tiles(indx1).sectionID;
    MP{counter}.pId= tid1;
    MP{counter}.p  = M{mix,1};
    
    MP{counter}.qz = L3(l2(mix)).tiles(indx2).sectionID;
    MP{counter}.qId= tid2;
    MP{counter}.q  = M{mix,2};
    counter = counter + 1;
end
js = pairs2json(MP); % generate json blob to be ingested into point-match database
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ingest js into point matches database
fn = ['temp_' num2str(randi(100000)) '_' num2str(lix) '.json'];
fid = fopen(fn, 'w');
fwrite(fid,js);
fclose(fid);
urlChar = sprintf('%s/owner/%s/matchCollection/%s/matches/', ...
    pm.server, pm.owner, pm.match_collection);
cmd = sprintf('curl -X PUT --connect-timeout 30 --header "Content-Type: application/json" --header "Accept: application/json" -d "@%s" "%s"',...
    fn, urlChar);
[a, resp]= evalc('system(cmd)');
delete(fn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% [8] solve whole system

%% [9] generate the new collection







































