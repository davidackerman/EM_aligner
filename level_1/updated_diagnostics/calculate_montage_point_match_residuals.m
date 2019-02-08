function [ output_struct ] = calculate_montage_point_match_residuals(...
    rc, point_matches, varargin)
%% Generate statistics about tile montage residuals
% Calculates tile residuals per section for renderer collection rc using
% point_matches point matches and varargin. varargin should contain either the
% unique_z to be analyzed, or zstart and zend required to obtain unique_z.
% varargin should contain 1-3 arguments:
%    1:                       The input should be unique_z
%    2:                       The input should be zstart and zend, or unique_z and options
%    3:                       The input should be zstart, zend, and options
% options fields and their defaults:
%    outlier_deviation_for_residuals : 10 Max point match residual for tile
%                                         before being considered an outlier
%    min_points      : 10     Minimum number of points for input to load_point_matches
%    output_data_per_tile : true  Output values and ratios for each tile
%    dir_scratch : /scratch/ackermand Scratch directory
%    filter_point_matches: true
%    verbose     : true               Output status

% Output:
%    output_struct:           Contains all tile-tile mean residuals, median
%                             of mean tile residuals (averaged over all tile-tile pairs), number of
%                             outliers of the median of mean residuals, outlier tile ids, number of
%                             unconnected tiles and the ids of unconnected tiles median and max of mean tile
%                             residuals.
% Author: Khaled Khairy, David Ackerman
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check the input arguments
options = [];
if isempty(varargin)
    error('Not enough input arguments: need either unique_z, or zstart and zend.');
elseif length(varargin)==1
    if isstruct(varargin{1})
        error('Not enough input arguments: need either unique_z, or zstart and zend.');
    else
        unique_z = varargin{1};
    end
elseif length(varargin)==2
    if isstruct(varargin{2})
        unique_z = varargin{1};
        options = varargin{2};
        zstart = unique_z(1);
        zend = unique_z(end);
    else
        zstart = varargin{1};
        zend = varargin{2};
        [unique_z, ~, ~, ~, ~] = get_section_ids(rc, zstart, zend+1);
        z_too_large = (unique_z>=zend+1);
        unique_z(z_too_large) = [];
    end
elseif length(varargin)==3
    if isstruct(varargin{3})
        zstart = varargin{1};
        zend = varargin{2};
        options = varargin{3};
        [unique_z, ~, ~, ~, ~] = get_section_ids(rc, zstart, zend+1);
        z_too_large = (unique_z>=zend+1);
        unique_z(z_too_large) = [];
    end
end

webopts = weboptions('Timeout', 60);

floor_unique_z = floor(unique_z);
unique_merged_z = unique(floor_unique_z);

new_dir_scratch=false;
if ~isfield(options, 'outlier_deviation_for_residuals'), options.outlier_deviation_for_residuals = 10; end
if ~isfield(options, 'min_points'), options.min_points = 10; end
if ~isfield(options,'output_data_per_tile'), options.output_data_per_tile = true; end
if ~isfield(options, 'dir_scratch')
    new_dir_scratch=true;
    options.dir_scratch = [pwd '/scratch_' num2str(randi(10000)) '_' datestr(datetime('now'),'yyyymmdd_HHMMSS')];
    warning('Will create temporary scratch directory %s which will be cleaned after', options.dir_scratch);
end
if ~isfield(options,'filter_point_matches'), options.filter_point_matches = true; end
if ~isfield(options, 'verbose'), options.verbose = true; end
if ~isfield(options, 'pmopts')
    options.pmopts.NumRandomSamplingsMethod = 'Desired confidence';
    options.pmopts.MaximumRandomSamples = 3000;
    options.pmopts.DesiredConfidence = 99.5;
    options.pmopts.PixelDistanceThreshold = 1;
end
if ~isfield(options, 'store_residual_outlier_pair_information'), options.store_residual_outlier_pair_information = false; end

dir_current = pwd;
dir_scratch = [options.dir_scratch '/temp_' num2str(randi(3000000))];
kk_mkdir(dir_scratch);
cd(dir_scratch);

num_el = length(unique_merged_z);
if options.output_data_per_tile, all_residuals_vector = cell(num_el,1); end

all_residuals_median = zeros(num_el,1);
all_residuals_mean = zeros(num_el,1);
all_residuals_max = zeros(num_el,1);
all_residuals_pair_max = zeros(num_el,1);
all_tile_ids_residuals_pair_max = cell(num_el,1);
all_pm_max = zeros(num_el,1);
all_tile_ids_pm_max = cell(num_el,1);
all_residuals_variance = zeros(num_el,1);
all_tile_ids = cell(num_el,1);
all_mapping_from_section_map_indices_to_current_indicies=cell(num_el,1);
all_residuals_outlier_count = zeros(num_el,1);
all_residuals_outlier_percent = zeros(num_el,1);
all_residuals_outlier_tile_indices = cell(num_el,1);
all_residuals_outlier_tile_ids = cell(num_el,1);
all_unconnected_count = zeros(num_el,1);
all_unconnected_tile_indices = cell(num_el,1);
all_unconnected_tile_ids = cell(num_el,1);
if options.store_residual_outlier_pair_information
    outlier_pm_max_pair_data = cell(num_el,1);
end
% Loop over all unique z and print out progress
if options.verbose
    fprintf('Montage Residuals Progress:');
    fprintf(['\n' repmat('.',1,num_el) '\n\n']);
end
[zu, sID, ~, ~, ns] = get_section_ids(rc, zstart, zend);
options.nbrs = 0;
parfor z_index = 1:numel(unique_merged_z)
    %% Determine point-matches and residuals for this section
    % First: load point-matches and section into "L" (point-matches are in L's pm struct field)
    matching_indices = find(floor_unique_z==unique_merged_z(z_index));
    valid_zs = unique_z(matching_indices);
%    section_information = [];
%     section_information.zu = zu(matching_indices);
%     section_information.sID = sID(matching_indices);
%     section_information.ns = ns(matching_indices);
% 
%     [L]  = ...
%         load_point_matches(valid_zs(1), valid_zs(end), rc, point_matches, 0, ...
%         options.min_points, 0,inf);
         [T, map_id, tIds, z_val, r, c] = load_all_transformations(rc, valid_zs, options.dir_scratch);
         [M, adj, W, np] = system_solve_helper_load_point_matches(valid_zs, options, point_matches, map_id, sID(matching_indices), size(T,1), r, c);
    
         %     if options.filter_point_matches
%         if isfield(options, 'pmopts')
%             L.pm = filter_pm(L.pm, options.pmopts);
%         else
%             L.pm = filter_pm(L.pm);
%         end
%     end
%     
    % Second: generate point-match residuals from L.pm by looping through
    % all point matches, transforming them and calculating the mean residual
    % for each tile pair
    tile_residuals = cell(numel(tIds),1);
    current_section_pair_max = -1;
    current_section_pm_max = -1;
    for point_match_index = 1:size(M,1)
        adjacent_tile_1 = adj(point_match_index,1);
        adjacent_tile_2 = adj(point_match_index,2);
        point_matches_tile_1 = M{point_match_index,1};
        point_matches_tile_2 = M{point_match_index,2};
        point_matches_tile_1 = [point_matches_tile_1 ones(size(point_matches_tile_1,1),1)]*[reshape(T(adjacent_tile_1,:),3,2), [0;0;1]];  % Apply transformation
        point_matches_tile_2 = [point_matches_tile_2 ones(size(point_matches_tile_2,1),1)]*[reshape(T(adjacent_tile_2,:),3,2), [0;0;1]];  % Apply transformation
        all_residuals = sqrt((point_matches_tile_1(:,1)-point_matches_tile_2(:,1)).*(point_matches_tile_1(:,1)-point_matches_tile_2(:,1))  + (point_matches_tile_1(:,2)-point_matches_tile_2(:,2)).* (point_matches_tile_1(:,2)-point_matches_tile_2(:,2)));       
        residual = mean(all_residuals); 
        [current_pm_max, current_pm_max_i] = max(all_residuals);
        tile_residuals{adjacent_tile_1} = [tile_residuals{adjacent_tile_1} residual];  % Aggregate residuals for adjacent tile 1
        tile_residuals{adjacent_tile_2} = [tile_residuals{adjacent_tile_2} residual];  % Aggregate residuals for adjacent tile 2
        if options.store_residual_outlier_pair_information % Store information about high residuals
           if current_pm_max >= options.outlier_deviation_for_residuals
               [sorted_tIds, sorted_order] = sort({tIds{adjacent_tile_1}, tIds{adjacent_tile_2}});
               if sorted_order(1)==1
                  pm_max_xyz = [point_matches_tile_1(current_pm_max_i,1), point_matches_tile_1(current_pm_max_i,2), unique_merged_z(z_index)]; 
               else
                  pm_max_xyz = [point_matches_tile_2(current_pm_max_i,1), point_matches_tile_2(current_pm_max_i,2), unique_merged_z(z_index)];
               end
               outlier_pm_max_pair_data{z_index}{end+1,1} = sprintf('%0.1f %f %s %s %f %f', pm_max_xyz(3), current_pm_max, sorted_tIds{1}, sorted_tIds{2}, pm_max_xyz(1), pm_max_xyz(2));
           end
        end
        if residual>current_section_pair_max
            current_section_pair_max=residual;
            all_tile_ids_residuals_pair_max{z_index} = {tIds(adjacent_tile_1), tIds(adjacent_tile_2)};
        end
        if current_pm_max>current_section_pm_max
            current_section_pm_max=current_pm_max;
            all_tile_ids_pm_max{z_index} = {tIds(adjacent_tile_1), tIds(adjacent_tile_2)};
        end
    end
    %Necessary so mapping order is the same between area/perimeter calculations
    %and residuals
    split_zs = unique_z(floor_unique_z == unique_merged_z(z_index));
    rc_data_for_mapping = [];
    for split_z_index = 1:numel(split_zs)
        urlChar = sprintf('%s/owner/%s/project/%s/stack/%s/z/%.1f/tile-specs', ...
            rc.baseURL, rc.owner, rc.project, rc.stack,split_zs(split_z_index) );
        rc_data_for_mapping = [rc_data_for_mapping; webread(urlChar, webopts)];
    end
    [~, mapping_from_section_map_indices_to_current_indicies] = ismember({rc_data_for_mapping(:).tileId}, tIds);
    if numel({rc_data_for_mapping(:).tileId})==numel(tIds) && sum(mapping_from_section_map_indices_to_current_indicies==0)==0 %Then all ids are covered
        all_mapping_from_section_map_indices_to_current_indicies{z_index} = mapping_from_section_map_indices_to_current_indicies;
    else
        error('Different tile ids');
    end
    
    all_residuals_pair_max(z_index) = current_section_pair_max;
    all_pm_max(z_index) = current_section_pm_max;
    %% Determine residual outliers
    % Separate connected and unconnected tiles
    %tile_ids = {L.tiles.renderer_id};
    [~, do_tiles_appear_in_adj] =ismember((1:numel(tIds)),unique(adj(:)));
    if numel(tIds)>1 % Make sure that there should be any connected tiles
        unconnected_tiles = (do_tiles_appear_in_adj == 0);
    else
        unconnected_tiles = []; %Then just a single tile
    end
    all_tile_ids{z_index} = tIds;
    all_unconnected_count(z_index) = sum(unconnected_tiles);
    all_unconnected_tile_ids{z_index} = tIds(unconnected_tiles);
    all_unconnected_tile_indices{z_index} = find(unconnected_tiles);
    if options.output_data_per_tile, all_residuals_vector{z_index} = tile_residuals; end  % Store tile residuals for this section
    % Calculate median of mean tile residuals, and outliers
    only_greater_than = true;
    
    max_residuals = cellfun(@max,tile_residuals,'UniformOutput',false);  % all tile residuals for section zu1(z_index)
    empties = cellfun('isempty',max_residuals);
    max_residuals(empties) = {NaN};
    max_residuals = cell2mat(max_residuals);
    [all_residuals_median(z_index), all_residuals_mean(z_index), all_residuals_max(z_index), all_residuals_variance(z_index), all_residuals_outlier_count(z_index), all_residuals_outlier_percent(z_index), all_residuals_outlier_tile_indices{z_index}, all_residuals_outlier_tile_ids{z_index}] = ...
        calculate_statistics_and_outliers(max_residuals, options.outlier_deviation_for_residuals, tIds, 'fixed_cutoff', only_greater_than);
    if options.verbose, fprintf('\b|\n'); end
end
% Create output struct
if options.output_data_per_tile
    output_struct.values = all_residuals_vector;
    output_struct.median_of_means = all_residuals_median;
    output_struct.mean_of_means = all_residuals_mean;
    output_struct.max_of_means = all_residuals_max;
    output_struct.variance_of_means = all_residuals_variance;
    output_struct.max_of_pairs = all_residuals_pair_max;
    output_struct.tile_ids_of_pairs = all_tile_ids_residuals_pair_max;
    output_struct.max_of_pm = all_pm_max;
    output_struct.tile_ids_of_pm = all_tile_ids_pm_max;
    output_struct.all_tile_ids = all_tile_ids;
    output_struct.all_mapping_from_section_map_indices_to_current_indicies=all_mapping_from_section_map_indices_to_current_indicies;
    output_struct.outlier_count = all_residuals_outlier_count;
    output_struct.outlier_percent = all_residuals_outlier_percent;
    output_struct.outlier_tile_indices =all_residuals_outlier_tile_indices;
    output_struct.outlier_tile_ids = all_residuals_outlier_tile_ids;
    output_struct.unconnected_count = all_unconnected_count;
    output_struct.unconnected_tile_indices =all_unconnected_tile_indices;
    output_struct.unconnected_tile_ids = all_unconnected_tile_ids;
else
    output_struct.median_of_means = all_residuals_median;
    output_struct.mean_of_means = all_residuals_mean;
    output_struct.max_of_means = all_residuals_max;
    output_struct.variance_of_means = all_residuals_variance;
    output_struct.max_of_pairs = all_residuals_pair_max;
    output_struct.tile_ids_of_pairs = all_tile_ids_residuals_pair_max;
    output_struct.max_of_pm = all_pm_max;
    output_struct.tile_ids_of_pm = all_tile_ids_pm_max;
    output_struct.outlier_count = all_residuals_outlier_count;
    output_struct.outlier_percent = all_residuals_outlier_percent;
    output_struct.outlier_tile_ids = all_residuals_outlier_tile_ids;
    output_struct.unconnected_count = all_unconnected_count;
    output_struct.unconnected_tile_ids = all_unconnected_tile_ids;
end
if options.store_residual_outlier_pair_information % Store information about high residuals
    output_struct.outlier_pm_max_pair_data = outlier_pm_max_pair_data;
end
cd(dir_current);
if new_dir_scratch
    system(sprintf('rm -rf %s', options.dir_scratch));
end
end

