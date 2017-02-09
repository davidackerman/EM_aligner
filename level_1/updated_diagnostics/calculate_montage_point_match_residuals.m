function [ output_struct ] = calculate_montage_point_match_residuals(rc, point_matches, options, unique_z)
all_tile_residuals_vector = cell(numel(unique_z),1);
all_residuals_median = zeros(numel(unique_z),1);
all_residuals_number_of_outliers = zeros(numel(unique_z),1);
all_residuals_outliers_tile_ids = cell(numel(unique_z),1);
all_number_of_unconnected_tiles = zeros(numel(unique_z),1);
all_unconnected_tile_ids = cell(numel(unique_z),1);

% Print out status and loop through all unique zs
fprintf('Montage Residuals Progress:');
fprintf(['\n' repmat('.',1,numel(unique_z)) '\n\n']);
parfor z_index = 1:length(unique_z)
    %% %%% determine point-matches, solution and residuals for this section
    % First: load point-matches and section into "L" (point-matches are in L's pm struct field)
    if (z_index + options.nbrs+1)>numel(unique_z)
        point_match_zlast = unique_z(end);
    else
        point_match_zlast = unique_z(z_index + options.nbrs);
    end
    [L]  = ...
        load_point_matches(unique_z(z_index), point_match_zlast, rc, point_matches, options.nbrs, ...
        options.min_points, 0);
    
    tile_residuals = cell(numel(L.tiles),1);
    % Second: generate point-match residuals from L.pm by transforming them and taking the sum of squared
    % residuals
    tile_ids = {L.tiles.renderer_id};
    for point_match_index = 1:size(L.pm.M,1)
        adjacent_tile_1 = L.pm.adj(point_match_index,1);
        adjacent_tile_2 = L.pm.adj(point_match_index,2);
        point_matches_tile_1 = L.pm.M{point_match_index,1};
        point_matches_tile_2 = L.pm.M{point_match_index,2};
        point_matches_tile_1 = [point_matches_tile_1 ones(size(point_matches_tile_1,1),1)]*L.tiles(adjacent_tile_1).tform.T;  % apply transformation
        point_matches_tile_2 = [point_matches_tile_2 ones(size(point_matches_tile_2,1),1)]*L.tiles(adjacent_tile_2).tform.T;  % apply transformation
        residual = mean(sqrt((point_matches_tile_1(:,1)-point_matches_tile_2(:,1)).*(point_matches_tile_1(:,1)-point_matches_tile_2(:,1))  + (point_matches_tile_1(:,2)-point_matches_tile_2(:,2)).* (point_matches_tile_1(:,2)-point_matches_tile_2(:,2))));    %%%% sum of squared residuals
        tile_residuals{adjacent_tile_1} = [tile_residuals{adjacent_tile_1} residual];  % aggregate residuals for tile a1
        tile_residuals{adjacent_tile_2} = [tile_residuals{adjacent_tile_2} residual];
    end
    all_tile_residuals_vector{z_index} = tile_residuals;  % store tile residuals for this section
    %% determine residual outliers
    % for each tile, calculate the median of means of point-match residuals
    [~, do_tiles_appear_in_adj] =ismember((1:numel(L.tiles)),unique(L.pm.adj(:)));
    unconnected_tiles = (do_tiles_appear_in_adj == 0);
    all_number_of_unconnected_tiles(z_index) = sum(unconnected_tiles);
    all_unconnected_tile_ids{z_index} = tile_ids(unconnected_tiles);
    tile_residuals(unconnected_tiles)=[];
    tile_ids(unconnected_tiles) = [];
    [all_residuals_median(z_index), all_residuals_outliers_tile_ids{z_index},all_residuals_number_of_outliers(z_index)] = calculate_median_and_outliers(cellfun(@mean,tile_residuals),options.nstd,tile_ids, 'std');
    fprintf('\b|\n');
end
output_struct.MontageResiduals.residuals = all_tile_residuals_vector;
output_struct.MontageResiduals.median_of_means = all_residuals_median;
output_struct.MontageResiduals.number_of_outliers = all_residuals_number_of_outliers;
output_struct.MontageResiduals.outlier_tile_ids = all_residuals_outliers_tile_ids;
output_struct.MontageResiduals.number_of_unconnected_tiles = all_number_of_unconnected_tiles;
output_struct.MontageResiduals.unconnected_tile_ids = all_unconnected_tile_ids;

end

