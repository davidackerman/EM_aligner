clc;
kk_clock;
% configure source collection
rcsource.stack          = 'v12_acquire_merged';
rcsource.owner          ='flyTEM';
rcsource.project        = 'FAFB00';
rcsource.service_host   = '10.37.5.60:8080';
rcsource.baseURL        = ['http://' rcsource.service_host '/render-ws/v1'];
rcsource.verbose        = 1;

% configure montage collection

rctarget_montage.stack          = ['EXP_v12_montage_22_03'];
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
opts.min_tiles = 2; % minimum number of tiles that constitute a cluster to be solved. Below this, no modification happens
opts.degree = 1;    % 1 = affine, 2 = second order polynomial, maximum is 3
opts.outlier_lambda = 1e3;  % large numbers result in fewer tiles excluded

opts.solver = 'backslash';
opts.min_points = 5;
opts.nbrs = 0;
opts.xs_weight = 1/100;
opts.stvec_flag = 0;   % i.e. do not assume rcsource providing the starting values.


% % test for best regularization parameter
% % This is the smallest that does not cause downscaling of tiles
% err = [];
% count = 1;
% for explambda = [-4: 5]
%     disp(explambda);
%     opts.lambda = 10^explambda;
%     opts.edge_lambda = 10^(explambda+1);
%     [mL2, A, err_res, R]= solve_slab(rcsource, pm, 22, 22, [], opts);
%     err(count) = err_res{1};
%     
%     %scl(count) = sum((mL2.tiles(20).tform.T([1 5])-1).^2)/2;
%     
%     scl(count) = sum(([mL2.tiles(20).tform.A(2) mL2.tiles(20).tform.A(2)]-1).^2)/2;
%     
%     count = count + 1;
% end
% explambda = [-4: 5];
% lambda = 10.^(explambda);
% plot(explambda, mat2gray(err), 'LineWidth', 1.0, 'Color', [1 0 0]);
% hold on;plot(explambda, mat2gray(scl), 'LineWidth', 1.0,'Color', [0 0 1]);

%

opts.lambda = 1e-3;
opts.edge_lambda = 1e-3;
% % solve montage for first section
% [mL1, A, err21]= solve_slab(rcsource, pm, 21, 21, [], opts);
% ingest_section_into_LOADING_collection(mL1, rctarget_montage, rcsource, pwd, 1);
% solve montage for second section
[mL2, A, err22, R]= solve_slab(rcsource, pm, 22, 22, [], opts);
ingest_section_into_LOADING_collection(mL2, rctarget_montage, rcsource, pwd, 1);
resp = set_renderer_stack_state_complete(rctarget_montage);

disp([ err22{1}]);


%% render
im = {};
% z = 21;
%  [Wbox, bbox, url] = get_section_bounds_renderer(rctarget_montage, z);
% [im{1}, v, url] = get_image_box_renderer(rctarget_montage, 21, Wbox, 0.1);
% figure;imshow(im{1});title(num2str(z));

z = 22;
 [Wbox, bbox, url] = get_section_bounds_renderer(rctarget_montage, z);
[im{2}, v, url] = get_image_box_renderer(rctarget_montage, z, Wbox, 0.1);
figure;imshow(im{2});title(num2str(z));