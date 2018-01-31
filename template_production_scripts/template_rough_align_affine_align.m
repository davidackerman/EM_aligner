%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% "stackId" : {
% 
%     "owner" : "flyTEM",
% 
%     "project" : "trautmane_fafb_fold_montage_tiles_sample_tier_0",
% 
%     "stack" : "0001x0001_000000"
% 
%   }
% 
%  
% 
% http://tem-services.int.janelia.org:8080/render-ws/view/stacks.html?dynamicRenderHost=renderer%3A8080&catmaidHost=renderer-catmaid%3A8000&renderStackOwner=flyTEM&matchOwner=flyTEM&renderStackProject=trautmane_fafb_fold_montage_tiles_sample_tier_0&renderStack=0001x0001_000000&matchCollection=trautmane_fafb_fold_montage_tiles_sample_tier_0_0001x0001_000000
% 
%  
% 
% The point match collection is:
% 
% "matchCollectionId" : {
% 
%       "owner" : "flyTEM",
% 
%       "name" : "trautmane_fafb_fold_montage_tiles_sample_tier_0_0001x0001_000000"
% 
%     }

% configure 
diary on;
clc
nfirst= 2200;%
nlast = 2699;
% configure source
% configure source collection
rcsource.stack          = '0001x0001_000000';
rcsource.owner          ='flyTEM';
rcsource.project        = 'trautmane_fafb_fold_montage_tiles_sample_tier_0';
rcsource.service_host   = '10.40.3.162:8080';
rcsource.baseURL        = ['http://' rcsource.service_host '/render-ws/v1'];
rcsource.verbose        = 0;


% configure rigid_approximation alignment
rcrigid.stack          = ['0001x0001_000000_rigid' num2str(nfirst) '_' num2str(nlast)];
rcrigid.owner          ='flyTEM';
rcrigid.project        = 'trautmane_fafb_fold_montage_tiles_sample_tier_0';
rcrigid.service_host   = '10.40.3.162:8080';
rcrigid.baseURL        = ['http://' rcrigid.service_host '/render-ws/v1'];
rcrigid.verbose        = 0;

% configure fine output collection
rcfine.stack          = ['0001x0001_000000_affine'];
rcfine.owner          ='flyTEM';
rcfine.project        = 'trautmane_fafb_fold_montage_tiles_sample_tier_0';
rcfine.service_host   = '10.40.3.162:8080';
rcfine.baseURL        = ['http://' rcfine.service_host '/render-ws/v1'];
rcfine.verbose        = 0;
rcfine.versionNotes   = 'Fine alignment: regularized by rigid';


% configure point-match collection
pm.server = 'http://10.40.3.162:8080/render-ws/v1';
pm.owner = 'flyTEM';
pm.match_collection = 'trautmane_fafb_fold_montage_tiles_sample_tier_0_0001x0001_000000';
pm.verbose = 0;

dir_scratch = '/groups/flyem/data/khairy_alignments/scratch/rigid_approximation_experiments_00';%'/scratch/khairyk';
kk_mkdir(dir_scratch);
cd(dir_scratch);
kk_clock();
% 

%%%  configure 
opts.degree = 0;    % 1 = affine, 2 = second order polynomial, maximum is 3
opts.solver = 'backslash';

opts.transfac = 1;  % translation parameter regidity
% opts.xfac = 1;   % 2nd order parameter rigidity in x
% opts.yfac = 1;   % 2nd order parameter regidity in y
% opts.lambda = 10^(6.5); % 10^4.5 (best results so far for affine) ------------------>

opts.nbrs = 2;
opts.nbrs_step = 1;
opts.xs_weight = 1.0;
opts.min_points = 1;
opts.max_points = 150;


opts.outlier_lambda = 1e2;  % large numbers result in fewer tiles excluded
opts.min_tiles = 2; % minimum number of tiles that constitute a cluster to be solved. Below this, no modification happens
opts.matrix_only = 0;   % 0 = solve , 1 = only generate the matrix
opts.distribute_A = 16;  % # shards of A
opts.dir_scratch = dir_scratch;

% opts.stvec_flag = 1;   % 0 = regularization against rigid model (i.e.; starting value is not supplied by rc)
opts.distributed = 0;
opts.disableValidation = 1;
%opts.edge_lambda = opts.lambda;
opts.use_peg = 0;

% % configure point-match filter
opts.filter_point_matches = 1;
opts.pmopts.NumRandomSamplingsMethod = 'Desired confidence';
opts.pmopts.MaximumRandomSamples = 5000;
opts.pmopts.DesiredConfidence = 99.9;
opts.pmopts.PixelDistanceThreshold = 50.0;
opts.verbose = 0;
opts.debug = 0;


%     
% % % produce a translation-only stack
disp('Solve for rigid');
rcrigid.versionNotes = gen_versionNotes(opts);
[err,R, Tout_translation, PM, Diagnostics] = system_solve_rigid_approximation(nfirst, nlast, rcsource, pm, opts, rcrigid);

% h = figure;plot([Diagnostics.rms_o(:,1) Diagnostics.rms(:,1)]);
% title('Comparison rms per tile (or section) between source and translation');
% xlabel('section number');ylabel('rms residual error');legend('Acquire', 'Rigid');
% disp('Done!');

%% configure Affine fine alignment
rcrough = rcrigid;
clear opts;
% configure solver
opts.transfac = 10^(-10); 
opts.lambda = 10^(6); % 

opts.constrain_by_z = 0;
opts.sandwich = 0;
opts.constraint_fac = 1e15;
opts.filter_point_matches = 1;
opts.save_matrix = 0;

% configure solver
opts.degree = 1;    % 1 = affine, 2 = second order polynomial, maximum is 3
opts.solver = 'backslash';%

opts.nbrs = 2;
opts.nbrs_step = 1;
opts.xs_weight = 1.0;
opts.min_points = 1;
opts.max_points = inf;


opts.outlier_lambda = 1e2;  % large numbers result in fewer tiles excluded
opts.min_tiles = 2; % minimum number of tiles that constitute a cluster to be solved. Below this, no modification happens
opts.pastix.ncpus = 8;
opts.pastix.parms_fn = '/nrs/flyTEM/khairy/FAFB00v13/matlab_production_scripts/params_file.txt';
opts.pastix.split = 1; % set to either 0 (no split) or 1
opts.matrix_only = 0;   % 0 = solve , 1 = only generate the matrix
opts.distribute_A = 16;  % # shards of A
opts.dir_scratch = '/scratch/khairyk';

% opts.stvec_flag = 1;   % 0 = regularization against rigid model (i.e.; starting value is not supplied by rc)
opts.distributed = 0;
opts.disableValidation = 1;
opts.edge_lambda = opts.lambda;
opts.use_peg = 0;

% % configure point-match filter
opts.filter_point_matches = 1;
opts.pmopts.NumRandomSamplingsMethod = 'Desired confidence';
opts.pmopts.MaximumRandomSamples = 5000;
opts.pmopts.DesiredConfidence = 99.9;
opts.pmopts.PixelDistanceThreshold = 50;
opts.verbose = 1;
opts.debug = 0;

% %
disp('---- Solve for affine -----');
rcfine.versionNotes = gen_versionNotes(opts);
[err,R, Tout_affine, D_affine] = system_solve_affine_with_constraint(nfirst, nlast, rcrough, pm, opts, rcfine);
disp(err);


%%

%%%%%%%%%%%%%%%%%%%%%%% uncomment for various diagnostics
% 
% 
%    disp('Aligned:');
%     disp(rcfine);
%     % % generate x y residual plots
% %    h = figure;plot([Diagnostics.rms_o(:,1) Diagnostics.rms(:,1) D_affine.rms(:,1)]);
%     h = figure;plot([ D_affine.rms(:,1)]);
% 
%     merr = median(D_affine.rms(:,1));
%     title(['PM with scale 0.4. Comparison rms per tile: acquire, translation and affine. median rms error affine: ' num2str(merr)]);
%     xlabel('tile number');ylabel('total pixel residual error');legend('Affine');
%    
%     figure;plot(D_affine.rms);
%     title(['PM with scale 0.4. median rms error affine: ' num2str(merr)]);
%     xlabel('section number');ylabel('rms in pixels');legend('Affine');
% 
%     % % % deformation-based diagnostics
%     close all;
%     diagnostics_opts.show_deformation = 0;
%     diagnostics_opts.show_residuals   = 0;
%     diagnostics_opts.residual_info    = 0;
%     diagnostics_opts.show_deformation_summary = 0;
%     [mA, mS, sctn_map, confidence, tile_areas, tile_perimeters, tidsvec,...
%     section_conf, residual_outliers_tid, area_outliers_tid, outliers_tid,  T] =...
%     gen_diagnostics(rcsource, rcfine, nfirst, nlast, [], diagnostics_opts);
% 
%     figure;hist(mA, 100);title('Histogram of tile-area ratio: deformed/original');

    %% generate cross section
%     [Wbox, bbox, url] = get_slab_bounds_renderer(rcfine);
% scale = 0.1;
% dx = 5/scale;
% x = Wbox(1) + Wbox(3)/2;
% y = Wbox(2);
% height = Wbox(4);
% disp([x y height dx]);
% [Iyz, Io] = get_yz_image_renderer(rcfine, x, y, dx, height, scale,...
%                                 nfirst, nlast, [2 2 1]);
% % imtool(Iyz); %   
% % zrange = nfirst:nlast;
% % I = imresize(Iyz,[height numel(zrange)*(zres/yres)]);
% 
% % %
% h =  size(Iyz,1);%20000;
% w =  size(Iyz,2);
% im = Iyz(1:h,:);
% I = imresize(im,[h/32 (w)]);
% figure(4);clf;imshow(I); title('Affine');
% drawnow;



% % %% compare to acquire
%     [Wbox, bbox, url] = get_slab_bounds_renderer(rcsource);
% scale = 0.1;
% dx = 5/scale;
% x = Wbox(1) + Wbox(3)/2;
% y = Wbox(2);
% height = Wbox(4);
% disp([x y height dx]);
% [Iyz, Io] = get_yz_image_renderer(rcsource, x, y, dx, height, scale,...
%                                 nfirst, nlast, [2 2 1]);
% % imtool(Iyz); %   
% % zrange = nfirst:nlast;
% % I = imresize(Iyz,[height numel(zrange)*(zres/yres)]);
% 
% % %
% h =  size(Iyz,1);%20000;
% w =  size(Iyz,2);
% im = Iyz(1:h,:);
% Iacq = imresize(im,[h/32 (w)]);
% figure(3);clf;imshow(Iacq); title('acquire-only');


%    [I, Io] = get_xz_image_renderer(rc, x, y, width, dy, scale, zstart, zfinish, res);
    %% generate mini stack (optional)
%     [Wbox, bbox, url, minZ, maxZ] = get_slab_bounds_renderer(rcfine);
%     
%     n_spark_nodes = 2;
%     bill_to = 'hess';
%     spark_dir = '/groups/flyem/data/render/spark_output';
%     dir_out_ministack = [dir_out '/' stack{ix}.name '_ministack'];% '/groups/flyTEM/home/khairyk/mwork/FIBSEM/mini_stacks'; % /groups/flyem/data/khairy_alignments/D08_09_ministack
%     max_images = 3000;
%     scale = 1.0;
%     zscale = 1.0;
%     minz = nfirst;
%     maxz = nlast;
%     minx = Wbox(1);%11405;
%     width = Wbox(3);% 2000
%     str = generate_mini_stack(rcfine, scale, zscale, dir_out, minz, maxz, minx, width,  n_spark_nodes, bill_to);
% 




%% 
% % % slopes and affine paraemters
% % rc.stack          = ['column_27_fine_pastix08']; % -------------->
% % rc.owner          ='hessh';
% % rc.project        = 'fly_em_pre_iso';
% % rc.service_host   = '10.40.3.162:8080';
% % rc.baseURL        = ['http://' rc.service_host '/render-ws/v1'];
% % rc.verbose        = 0;
% % 
% % pm.server = 'http://10.40.3.162:8080/render-ws/v1';
% % pm.owner = 'hessh';
% % pm.match_collection = 'fly_em_pre_iso_column_27';
% % pm.verbose = 0;
% % 
% % zs = get_section_ids(rc);
% % 
% % slopes = zeros(numel(zs),4);
% % transforms = zeros(numel(zs),9);
% % fprintf(['\n' repmat('.',1,numel(zs)) '\n\n']);
% % parfor i=1:numel(zs)
% %     [x1, y1, ~, L1, ~]=get_tile_centers(rc, zs(i));
% %     allSlopes = (y1(2:end)-y1(1:end-1))./(x1(2:end)-x1(1:end-1));
% %     if numel(y1)>1
% %     slopes(i,:) = [ mean((y1(2:end)-y1(1:end-1))./(x1(2:end)-x1(1:end-1))), max(allSlopes), min(allSlopes), (y1(end)-y1(1))/(x1(end)-x1(1))];
% %     else
% %            slopes(i,:) = [ NaN, NaN, NaN, NaN]; 
% %     end
% %     totalTransforms = zeros(3);
% %     for j=1:numel(L1.tiles)
% %         totalTransforms = totalTransforms + L1.tiles(j).tform.T;
% %     end
% %     totalTransforms = totalTransforms/numel(L1.tiles);
% %     transforms(i,:) = totalTransforms(:);
% %     fprintf('\b|\n');
% % end
% % figure();
% % plot(zs,slopes)
% % ylim([-0.004,0.004])
% % legend({'mean slope','max slope', 'min slope', 'start to end slope'})
% % 
% % ylabel('Slope');
% % xlabel('Z Value');
% % figure();
% % orderToPlot = [1, 4, 7,...
% %                2, 5, 8,...
% %                3, 6, 9];
% % for i=1:9
% %    subplot(3,3,i);
% %    plot(zs,transforms(:,orderToPlot(i)));
% %    title(['T(' num2str(ceil(i/3)) ',' num2str(mod(i-1,3)+1) ')']);
% %    grid on;
% % end
% % 
% % % %% montages
% % % opts.min_points = 1;
% % % opts.max_points = Inf;
% % % all_outputs = [];
% % % for i=24:27
% % %     if i==24
% % %         rc.stack = 'column_24_fine_03';
% % %     elseif i==25
% % %         rc.stack = 'column_25_fine_07';
% % %     elseif i==26
% % %         rc.stack = 'column_26_fine_pastix05';
% % %     else
% % %         rc.stack = 'column_27_fine_pastix06';
% % %     end
% % %         zs = get_section_ids(rc);
% % %     pm.match_collection = ['fly_em_pre_iso_column_' num2str(i)];
% % %     output = calculate_montage_point_match_residuals(rc, pm, zs(1), zs(end), opts);
% % %     output.zs = zs;
% % %     if i==24
% % %         all_outputs = output;
% % %     else
% % %         all_outputs(i-23) = output;
% % %     end
% % % end
