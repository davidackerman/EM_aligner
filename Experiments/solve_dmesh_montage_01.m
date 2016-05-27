clc;
kk_clock;


n = 22;

% % solve for rigid and store in database
% configure source collection
rcsource.stack          = 'v12_acquire_merged';
rcsource.owner          ='flyTEM';
rcsource.project        = 'FAFB00';
rcsource.service_host   = '10.37.5.60:8080';
rcsource.baseURL        = ['http://' rcsource.service_host '/render-ws/v1'];
rcsource.verbose        = 1;

% configure montage collection

rctarget_montage.stack          = ['EXP_v12_montage_P0_L1'];
rctarget_montage.owner          ='flyTEM';
rctarget_montage.project        = 'test';
rctarget_montage.service_host   = '10.37.5.60:8080';
rctarget_montage.baseURL        = ['http://' rctarget_montage.service_host '/render-ws/v1'];
rctarget_montage.verbose        = 1;

% configure point-match collection
pm.server           = 'http://10.40.3.162:8080/render-ws/v1';
pm.owner            = 'flyTEM';
pm.match_collection = 'v12_dmesh';

% configure solver
opts.min_tiles = 10; % minimum number of tiles that constitute a cluster to be solved. Below this, no modification happens
opts.degree = 0;    % 1 = affine, 2 = second order polynomial, maximum is 3
opts.outlier_lambda = 1e3;  % large numbers result in fewer tiles excluded
opts.solver = 'backslash';
opts.min_points = 5;
opts.nbrs = 0;
opts.xs_weight = 1/100;
opts.stvec_flag = 0;   % i.e. do not assume rcsource providing the starting values.
opts.distibuted = 0;

%
opts.lambda = 1e-3;
opts.edge_lambda = 1e-3;
% solve montage
[mL2, A, err22, R]= solve_slab(rcsource, pm, n, n, [], opts);
ingest_section_into_LOADING_collection(mL2, rctarget_montage, rcsource, pwd, 1);
resp = set_renderer_stack_state_complete(rctarget_montage);

%% fast solve for affine based on rigid collection
% configure source collection
rcsource.stack          = 'EXP_v12_montage_P0_L1';
rcsource.owner          ='flyTEM';
rcsource.project        = 'test';
rcsource.service_host   = '10.37.5.60:8080';
rcsource.baseURL        = ['http://' rcsource.service_host '/render-ws/v1'];
rcsource.verbose        = 1;

% configure point-match collection
pm.server           = 'http://10.40.3.162:8080/render-ws/v1';
pm.owner            = 'flyTEM';
pm.match_collection = 'v12_dmesh';

% configure solver
opts.min_tiles = 10; % minimum number of tiles that constitute a cluster to be solved. Below this, no modification happens
opts.degree = 1;    % 1 = affine, 2 = second order polynomial, maximum is 3
opts.outlier_lambda = 1e3;  % large numbers result in fewer tiles excluded
opts.solver = 'backslash';
opts.min_points = 5;
opts.nbrs = 0;
opts.xs_weight = 1/100;
opts.stvec_flag = 1;   % 0 =  do not assume that  rcsource provides the starting values.
opts.distibuted = 0;

%
opts.lambda = 1e-3;
opts.edge_lambda = 1e-3;
tdim = 6;
[K, Lm, tId, T, To, sectionId, z] = matrix_system_gen(n, n, ...
                                            rcsource, pm, opts.nbrs, ...
                                            opts.min_points, opts.xs_weight,...
                                            opts.stvec_flag, tdim, opts.lambda);

                                        
% read json tile specs
urlChar = sprintf('%s/owner/%s/project/%s/stack/%s/z/%d/last-tile-transforms', ...
    rcsource.baseURL, rcsource.owner, rcsource.project, rcsource.stack, n);
jt = webread(urlChar);
% read point-matches
urlChar = sprintf('%s/owner/%s/matchCollection/%s/group/%s/matchesWithinGroup', ...
            pm.server, pm.owner, pm.match_collection, '1.0');
jp = webread(urlChar);

































