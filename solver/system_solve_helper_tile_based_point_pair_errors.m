function [tile_err] = system_solve_helper_tile_based_point_pair_errors(pm, res, ntiles)

tpr = cell(ntiles,1);

% initialize
for ix = 1:ntiles
    tpr{ix} = [];
end


% aggregate absolute point-pair error for each tile
count = 1;
for ix = 1:size(pm.M,1)
    indx1 = pm.adj(ix,1);
    indx2 = pm.adj(ix,2);
    
    vecu = count:count+pm.np(ix)-1;
    vecv = count+pm.np(ix):count+2*pm.np(ix)-1;
    npoints = numel(vecu);
    rx = sum(abs(res(vecu)),1)/npoints;
    ry = sum(abs(res(vecv)),1)/npoints;
    tpr{indx1} = [tpr{indx1};rx ry];
    tpr{indx2} = [tpr{indx2};rx ry];
    dr = sqrt(rx.*rx + ry.*ry);

    
    count = count + 2*pm.np(ix);
end

% average the error
delix = zeros(ntiles,1, 'logical');
tile_err = zeros(ntiles,2);
parfor ix = 1:ntiles
    if isempty(tpr{ix})
        tile_err(ix,:) = [-999 -999];%full(sum(tpr{ix},1)/size(tpr{ix},1));
        delix(ix) = 1;
    else
        tile_err(ix,:) = full(sum(tpr{ix},1)/size(tpr{ix},1));
    end
end
tile_err(delix,:) = [];
