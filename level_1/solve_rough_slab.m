function [L2, needs_correction, pmfn, zsetd, zrange, t,dir_spark_work, cmd_str, fn_ids, ...
    target_solver_path, target_ids, target_matches, target_layer_images] = ...
    ...
    solve_rough_slab(dir_store_rough_slab, rcsource, rctarget_montage, rctarget_rough,ms, nfirst, nlast, dir_rough_intermediate_store, ...
    run_now, precalc_ids, precalc_matches, precalc_path)
%
% Author: Khaled Khairy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [0] configure collections and prepare quantities
% clc; clear all;
% kk_clock;
% nfirst = 1;
% nlast  = 15;
%
% % configure source collection
% rcsource.stack          = 'v12_acquire_merged';
% rcsource.owner          ='flyTEM';
% rcsource.project        = 'FAFB00';
% rcsource.service_host   = '10.37.5.60:8080';
% rcsource.baseURL        = ['http://' rcsource.service_host '/render-ws/v1'];
% rcsource.verbose        = 1;
%
% % configure montage collection
% rctarget_montage.stack          = ['EXP_dmesh_montage_P1'];
% rctarget_montage.owner          ='flyTEM';
% rctarget_montage.project        = 'test';
% rctarget_montage.service_host   = '10.37.5.60:8080';
% rctarget_montage.baseURL        = ['http://' rctarget_montage.service_host '/render-ws/v1'];
% rctarget_montage.verbose        = 1;
%
%
% % configure rough collection
% rctarget_rough.stack          = ['EXP_dmesh_rough_P1_' num2str(nfirst) '_' num2str(nlast)];
% rctarget_rough.owner          ='flyTEM';
% rctarget_rough.project        = 'test';
% rctarget_rough.service_host   = '10.37.5.60:8080';
% rctarget_rough.baseURL        = ['http://' rctarget_rough.service_host '/render-ws/v1'];
% rctarget_rough.verbose        = 1;
%
%
% % configure montage-scape point-match generation
% ms.service_host                 = rctarget_montage.service_host;
% ms.owner                        = rctarget_montage.owner;
% ms.project                      = rctarget_montage.project;
% ms.stack                        = rctarget_montage.stack;
% ms.first                        = num2str(nfirst);
% ms.last                         = num2str(nlast);
% ms.fd_size                      = '8';
% ms.min_sift_scale               = '0.55';
% ms.max_sift_scale               = '1.0';
% ms.steps                        = '3';
% ms.scale                        = '0.1';    % normally less than 0.05 -- can be large (e.g. 0.2) for very small sections (<100 tiles)
% ms.similarity_range             = '15';
% ms.skip_similarity_matrix       = 'y';
% ms.skip_aligned_image_generation= 'y';
% ms.base_output_dir              = '/nobackup/flyTEM/spark_montage';
% ms.run_dir                      = ['scale_' ms.scale];
% ms.script                       = '/groups/flyTEM/home/khairyk/EM_aligner/renderer_api/generate_montage_scape_point_matches.sh';%'../unit_tests/generate_montage_scape_point_matches_stub.sh'; %
% ms.number_of_spark_nodes        = '2.0';

[zu, sID, sectionId, z, ns] = get_section_ids(rcsource, nfirst, nlast);

%% define target directores and file names for storage
%dir_rough_intermediate_store = '/nobackup/flyTEM/khairy/FAFB00v13/montage_scape_pms';
diary on;
target_solver_path = [dir_rough_intermediate_store '/solver_' num2str(nfirst) '_' num2str(nlast) ];
target_ids = [target_solver_path '/ids.txt'];
target_matches = [target_solver_path '/matches.txt'];
target_layer_images = [target_solver_path];

if run_now==0
    precalc_ids = target_ids;
    precalc_matches = target_matches;
    precalc_path = target_solver_path;
end
    
    
%% [2] START SPARK: generate montage-scapes and montage-scape point-matches
if run_now
    [L2, needs_correction, pmfn, zsetd, zrange, t,dir_spark_work,...
        cmd_str, fn_ids, missing_images, existing_images] = ...
        generate_montage_scapes_SIFT_point_matches(ms, run_now);
else
    [L2, needs_correction, pmfn, zsetd, zrange, t,dir_spark_work,...
        cmd_str, fn_ids, missing_images, existing_images] = ...
        generate_montage_scapes_SIFT_point_matches(ms, run_now, precalc_ids, precalc_matches, precalc_path);
end

%% organize directories/files and store intermediate data
source_layer_images = [dir_spark_work '/layer_images'];
if run_now
kk_mkdir(target_layer_images);
movefile(source_layer_images, target_layer_images);
movefile(pmfn, target_matches);
movefile(fn_ids, target_ids);
end
[zu, sID, sectionId, z, ns, zuf] = get_section_ids(rcsource, nfirst, nlast);
zuf(find(intersect(zu, missing_images))) = [];

for imix = 1:numel(zuf)%str2double(ms.first):str2double(ms.last)
    fn_im = sprintf('%s/layer_images/%.1f.png', target_solver_path,zuf(imix));
    L2.tiles(imix).path = fn_im;
    t(imix).path = fn_im;
%     counter = counter + 1;
end

if needs_correction==0
    
    %% [2'] correct if needed
    % if needs_correction   %%%%%%%%%% manually correct
    %     disp('Please interrup <ctrl-c> and proceed to manually set point correspondences');
    %     disp('generate images and then use cpselect(im1,im2)');
    %     disp('When done, set run-now to zeros and re-run generate_montage_scapes_SIFT_point_matches');
    %     pause;
    %
    %     [~, ai] = sort(zrange(:,2),'ascend');
    %     disp('Available connected tiles:');
    %     disp(zrange(ai,:));
    %
    %     % define point-match pairs that need to be connected
    %     c = [1 2; 18 19;26 27];
    %
    %     p = 3;
    %     im1 = imread(t(c(p,1)).path);
    %     im2 = imread(t(c(p,2)).path);
    %     cpselect(im1,im2);
    %
    %     im1p = movingPoints4;
    %     im2p = fixedPoints4;
    %     %%% try SURF first
    % %         p1  = detectSURFFeatures(im1, 'NumOctaves', 10,...
    % %                                    'NumScaleLevels', 5,...
    % %                                    'MetricThreshold', 1000);
    % %
    % %     [f1, vp1]  = extractFeatures(im1,  p1);
    % %     [f2, vp2]  = extractFeatures(im2,  p1);
    % %     [m1, m2, ~]  = im_pair_match_features(f1, vp1, f2, vp2);
    % %     showMatchedFeatures(im1,im2,im1p,im2p, 'montage', 'PlotOptions', {'ro', 'g+', 'w-'});h2 = findobj('Type', 'line');for lix = 1:numel(h2), set(h2(lix), 'LineWidth',2);end
    % %
    %
    %     %%% use cpselect for each image pair manually and append point-match file with new point-matches
    %         fid = fopen(pmfn, 'a+');
    %
    %      str = '';
    %     for pix = 1:size(im1p,1)
    %     str = [str sprintf('%.0f\t%.12f\t%.12f\t%.0f\t%.12f\t%.12f\n', ...
    %         t(c(p,1)).id, im1p(pix,1),im1p(pix,2), t(c(p,2)).id, im2p(pix,1),im2p(pix,2))];
    %     end
    %
    %     fprintf(fid, '%s', str);
    %     fclose(fid);
    %
    %
    %
    %     % manually run all corrections then re-run point-match generation
    %     run_now = 0;
    %     [L2, needs_correction, pmfn, zsetd, zrange, t,dir_spark_work, cmd_str,...
    %         fn_ids] = generate_montage_scapes_SIFT_point_matches(ms, run_now);
    %
    %
    % end
    %%
    % disp('Filtering point-matches using RANSAC');
    %
    % % %% filter point matches using RANSAC
    % geoTransformEst = vision.GeometricTransformEstimator; % defaults to RANSAC
    % geoTransformEst.Method = 'Random Sample Consensus (RANSAC)';%'Least Median of Squares';
    % geoTransformEst.Transform = 'Nonreflective similarity';%'Affine';%
    % geoTransformEst.NumRandomSamplingsMethod = 'Desired confidence';
    % geoTransformEst.MaximumRandomSamples = 1500;
    % geoTransformEst.DesiredConfidence = 99.99;
    % for pmix = 1:size(L2.pm.M,1)
    %     m1 = L2.pm.M{pmix,1};
    %     m2 = L2.pm.M{pmix,2};
    %     % Invoke the step() method on the geoTransformEst object to compute the
    %     % transformation from the |distorted| to the |original| image. You
    %     % may see varying results of the transformation matrix computation because
    %     % of the random sampling employed by RANSAC.
    %     [tform_matrix, inlierIdx] = step(geoTransformEst, m2, m1);
    %     m1 = m1(inlierIdx,:);
    %     m2 = m2(inlierIdx,:);
    %     L2.pm.M{pmix,1} = m1;
    %     L2.pm.M{pmix,2} = m2;
    %     w = L2.pm.W{pmix};
    %     L2.pm.W{pmix} = w(inlierIdx);
    % end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    % % check point match quality
    % for lix = 1%:size(L2.pm.adj,1)
    %     ix1 = L2.pm.adj(lix,1);
    %     ix2 = L2.pm.adj(lix,2);
    %     t1 = L2.tiles(ix1);
    %     t2 = L2.tiles(ix2);
    %     M = L2.pm.M(lix,:);
    %     show_feature_point_correspondence(t1,t2,M);title([num2str(ix1) '   ' num2str(ix2)]);
    %     drawnow;
    % end
    %% [3] rough alignment solve for montage-scapes
    disp('Solving');
    % solve
    rigid_opts.apply_scaling = 1;
    [mLR, errR, mLS] = get_rigid_approximation(L2, 'backslash', rigid_opts);  % generate rigid approximation to use as regularizer
    mL3 = mLR;
%     affine_opts.lambda = 1e-3;
%     [mL3, errA] = solve_affine_explicit_region(mLR, affine_opts); % obtain an affine solution
%     
    
    mL3s = split_z(mL3);
    %% [4] apply rough alignment to montaged sections (L_montage) and generate "rough_aligned" collection  %% %%%%%% sosi
    disp('Apply rough alignment to full set of montaged sections:');
    indx = find(zu-floor(zu)>0);
    zu(indx) = [];
    diary off;
    diary on;
    % load montages
    disp('-- Loading montages');
    tic
    parfor zix = 1:numel(zu),
        L_montage(zix) = Msection(rctarget_montage, zu(zix));
    end
    toc
    %%%% apply rough alignment solution to montages
    disp('-- Retrieve bounding box for each section');
    tic
    parfor lix = 1:numel(L_montage), L_montage(lix) = get_bounding_box(L_montage(lix));end
    mL3.update_tile_info_switch = -1;
    mL3 = get_bounding_box(mL3);
    Wbox = [mL3.box(1) mL3.box(3) mL3.box(2)-mL3.box(1) mL3.box(4)-mL3.box(3)];disp(Wbox);
    wb1 = Wbox(1);
    wb2 = Wbox(2);
    toc
    disp('-- Perform transformation for all tiles in each section');
    tic
    fac = str2double(ms.scale); %0.25;
    smx = [fac 0 0; 0 fac 0; 0 0 1]; %scale matrix
    invsmx = [1/fac 0 0; 0 1/fac 0; 0 0 1];
    tmx2 = [1 0 0; 0 1 0; -wb1 -wb2 1]; % translation matrix for montage_scape stack
    for lix = 1:numel(L_montage)
         b1 = L_montage(lix).box;
         dx = b1(1);dy = b1(3);
        tmx1 = [1 0 0; 0 1 0; -dx -dy 1];  % translation matrix for section box
        for tix = 1:numel(L_montage(lix).tiles)
            newT = L_montage(lix).tiles(tix).tform.T * tmx1 * smx * mL3s(lix).tiles(1).tform.T * tmx2 * (invsmx);
            %newT = L3(lix).tiles(tix).tform.T * tmx1 * invsmx * mL3s(lix).tiles(1).tform.T * tmx1 * (smx);
            L_montage(lix).tiles(tix).tform.T = newT;
        end
        L_montage(lix) = update_XY(L_montage(lix));
    end
    opts.outlier_lambda = 1e3;  % large numbers result in fewer tiles excluded
    L_montage = concatenate_tiles(L_montage);
    toc
    
    % save 
    try
    fn = sprintf('%s//rough_aligned_slab_%d_%d.mat', dir_store_rough_slab, nfirst, nlast);
    save(fn, 'L_montage', 'rctarget_rough', 'rcsource');
    catch err_save
        kk_err_disp(err_save);
    end
    %%% ingest into renderer database as rough collection
        disp('-- Ingest into renderer database using overwrite.');
        diary off;
    diary on;

    tic
    ingest_section_into_renderer_database_overwrite(L_montage,rctarget_rough, rcsource, pwd);
    toc
    disp('Done!');
    diary off;
    
    
    
end



%      -0.98988    0.0097416            0
%     -0.011215     -0.99088            0
%         24068        10731            1
% 
%       -1.0082    0.0061571            0
%    -0.0056041      -1.0051            0
%         22461        10581            1
% 
%      -0.99381   -0.0044445            0
%     0.0067667     -0.99783            0
%         23843        10013            1
% 
% 
























