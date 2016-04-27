function [L, L_vec, pm_mx, err, scl, h] = solver_regularization_parameter_sweep(nfirst, nlast, rcsource, pm, opts, regstart, regfinish, step)
% Performs a regularization parameter sweep between regstart and regfinish (inclusive)
% Test for best regularization parameter graphically -- sosi: autosuggest in future
% This is the smallest parameter value that does not cause downscaling of tiles
% Author: Khaled Khairy. Janelia Research Campus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[L, ~, ~, pm_mx] = load_point_matches(nfirst, nlast, rcsource, pm, opts.nbrs, opts.min_points, opts.xs_weight); % disp(pm_mx{ix});
[L_vec, ntiles] = reduce_to_connected_components(L);
L_vec(ntiles<10) = [];

solve_time = [];
err = [];
scl = [];
count = 1;
for explambda = [regstart:step: regfinish]
    disp(explambda);
    kkt = tic;
    opts.lambda = 10^explambda;
    opts.edge_lambda = 10^(explambda);
    %[mL2, A, err_res, R]= solve_slab(rcsource, pm, nfirst, nlast, [], opts);
    [mL2, err_res, R] = solve_clusters(L_vec, opts, opts.stvec_flag);   % solves individual clusters and reassembles them into one
    err(count) = max(err_res);
    
    % measure deformation
    parfor ix = 1:numel(mL2.tiles)
        t = mL2.tiles(ix);
        [U S V] = svd(t.tform.T(1:2, 1:2));
        detS(ix) = det(S);
    end
%     states = [mL2.tiles(:).state]==1;
%     detS = detS(states); % remove entries for discarded tiles
    scl(count) = sum((detS-1).^2);  % for affine
    
    %scl(count) = sum(([mL2.tiles(20).tform.A(2) mL2.tiles(20).tform.A(2)]-1).^2)/2;  % for higher order
    %polynomials
    solve_time(count) = toc(kkt);
    disp(['time: ' num2str(solve_time(count))]);
    count = count + 1;
end

% plot
explambda = [regstart:step: regfinish];
h = figure;
plot(explambda, mat2gray(err), 'LineWidth', 1.0, 'Color', [1 0 0]);
hold on;
plot(explambda, mat2gray(scl), 'LineWidth', 1.0,'Color', [0 0 1]);