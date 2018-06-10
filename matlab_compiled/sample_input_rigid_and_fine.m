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
clear all
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
rcrigid.stack          = ['davidrigid_' num2str(nfirst) '_' num2str(nlast)];
rcrigid.owner          ='flyTEM';
rcrigid.project        = 'trautmane_fafb_fold_montage_tiles_sample_tier_0';
rcrigid.service_host   = '10.40.3.162:8080';
rcrigid.baseURL        = ['http://' rcrigid.service_host '/render-ws/v1'];
rcrigid.verbose        = 0;
 
% configure fine output collection
rcfine.stack          = ['davidTest_affine'];
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
 
%dir_scratch = '/groups/flyem/data/khairy_alignments/scratch/rigid_approximation_experiments_00';%'/scratch/khairyk';

%%%  configure 
rigid_solver_options.degree = 0;    % 1 = affine, 2 = second order polynomial, maximum is 3
rigid_solver_options.solver = 'backslash';
 
rigid_solver_options.transfac = 1;  % translation parameter regidity
% opts.xfac = 1;   % 2nd order parameter rigidity in x
% opts.yfac = 1;   % 2nd order parameter regidity in y
% opts.lambda = 10^(6.5); % 10^4.5 (best results so far for affine) ------------------>
 
rigid_solver_options.nbrs = 2;
rigid_solver_options.nbrs_step = 1;
rigid_solver_options.xs_weight = 1.0;
rigid_solver_options.min_points = 1;
rigid_solver_options.max_points = 150;
 
 
rigid_solver_options.outlier_lambda = 1e2;  % large numbers result in fewer tiles excluded
rigid_solver_options.min_tiles = 2; % minimum number of tiles that constitute a cluster to be solved. Below this, no modification happens
rigid_solver_options.matrix_only = 0;   % 0 = solve , 1 = only generate the matrix
rigid_solver_options.distribute_A = 16;  % # shards of A
rigid_solver_options.dir_scratch = '/scratch/ackermand';
 
% rigid_solver_options.stvec_flag = 1;   % 0 = regularization against rigid model (i.e.; starting value is not supplied by rc)
rigid_solver_options.distributed = 0;
rigid_solver_options.disableValidation = 1;
%rigid_solver_options.edge_lambda = rigid_solver_options.lambda;
rigid_solver_options.use_peg = 0;
 
% % configure point-match filter
rigid_solver_options.filter_point_matches = 1;
rigid_solver_options.pmopts.NumRandomSamplingsMethod = 'Desired confidence';
rigid_solver_options.pmopts.MaximumRandomSamples = 5000;
rigid_solver_options.pmopts.DesiredConfidence = 99.9;
rigid_solver_options.pmopts.PixelDistanceThreshold = 50.0;
rigid_solver_options.verbose = 0;
rigid_solver_options.debug = 0;

%% configure Affine fine alignment
% configure solver
fine_solver_options.transfac = 10^(-10); 
fine_solver_options.lambda = 10^(6); % 
 
fine_solver_options.constrain_by_z = 0;
fine_solver_options.sandwich = 0;
fine_solver_options.constraint_fac = 1e15;
fine_solver_options.filter_point_matches = 1;
fine_solver_options.save_matrix = 0;
 
% configure solver
fine_solver_options.degree = 1;    % 1 = affine, 2 = second order polynomial, maximum is 3
fine_solver_options.solver = 'backslash';%
 
fine_solver_options.nbrs = 2;
fine_solver_options.nbrs_step = 1;
fine_solver_options.xs_weight = 1.0;
fine_solver_options.min_points = 1;
fine_solver_options.max_points = inf;
 
 
fine_solver_options.outlier_lambda = 1e2;  % large numbers result in fewer tiles excluded
fine_solver_options.min_tiles = 2; % minimum number of tiles that constitute a cluster to be solved. Below this, no modification happens
fine_solver_options.pastix.ncpus = 8;
fine_solver_options.pastix.parms_fn = '/nrs/flyTEM/khairy/FAFB00v13/matlab_production_scripts/params_file.txt';
fine_solver_options.pastix.split = 1; % set to either 0 (no split) or 1
fine_solver_options.matrix_only = 0;   % 0 = solve , 1 = only generate the matrix
fine_solver_options.distribute_A = 16;  % # shards of A
fine_solver_options.dir_scratch = '/scratch/ackermand';
 
% opts.stvec_flag = 1;   % 0 = regularization against rigid model (i.e.; starting value is not supplied by rc)
fine_solver_options.distributed = 0;
fine_solver_options.disableValidation = 1;
fine_solver_options.edge_lambda = fine_solver_options.lambda;
fine_solver_options.use_peg = 0;
 
% % configure point-match filter
fine_solver_options.filter_point_matches = 1;
fine_solver_options.pmopts.NumRandomSamplingsMethod = 'Desired confidence';
fine_solver_options.pmopts.MaximumRandomSamples = 5000;
fine_solver_options.pmopts.DesiredConfidence = 99.9;
fine_solver_options.pmopts.PixelDistanceThreshold = 50;
fine_solver_options.verbose = 1;
fine_solver_options.debug = 0;
%%
sl.first_section = nfirst;
sl.last_section = nlast;
sl.rigid_solver_options = rigid_solver_options;
sl.fine_solver_options = fine_solver_options;
sl.source_collection = rcsource;
sl.rigid_collection = rcrigid;
sl.fine_collection = rcfine;
sl.source_point_match_collection = pm;
sl.verbose = 1;

fn = [pwd '/sample_rigid_and_fine.json'];
str = savejson('', sl);
fid = fopen(fn,'w');
fprintf(fid,str);
fclose(fid);
%% make sure we can read this file
options = loadjson(fileread(fn));
 