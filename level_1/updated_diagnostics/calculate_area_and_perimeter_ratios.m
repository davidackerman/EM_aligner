function [output_struct, varargout] = calculate_area_and_perimeter_ratios( rcsource, rc, zstart, zend, options )
%% generate statistics about residuals and tile deformation
% Summarizes point-match residuals and tile deformation per tile and section taking
% into accounts its neighbors.
% opts fields and their defaults:
%    min_points     : 5
%    nbrs           : 4
%    show_deformation: 1      0 = don'e show
%                             1 = display visible figure
%                             2 = save image of invisible figure to disk (not implemented yet)
%                             3 = captures invisible figure to memory (not implemented yet)
%    show_residuals: 1        0 = don't show
%                             1 = displays a visible figure
%                             2 = save image of invisiblle figure to disk (not implemented yet)
%                             3 = captures invisible figure to memory (not implemented yet)
% Output:
%       mA, mS
% Author: Khaled Khairy, David Ackerman
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<5
    options.nstd = 2;
end

[unique_z, section_ids_grouped_by_z, ~, ~, ~] = get_section_ids(rc, zstart, zend);

%% find corresponding tile centers
% initialize variables to store deformation and residuals
all_rc_section_map  = cell(numel(unique_z),1);
all_rc_area_ratio_outliers_tile_ids = cell(numel(unique_z),1);
all_rc_perimeter_ratio_outliers_tile_ids = cell(numel(unique_z),1);
all_rc_tile_areas = cell(numel(unique_z),1);
all_rc_tile_perimeters = cell(numel(unique_z),1);
all_rc_tile_area_ratios = cell(numel(unique_z),1);
all_rc_tile_perimeter_ratios = cell(numel(unique_z),1);
all_rc_heights = cell(numel(unique_z),1);
all_rc_widths = cell(numel(unique_z),1);
all_rc_area_ratio_median = zeros(numel(unique_z),1);
all_rc_perimeter_ratio_median = zeros(numel(unique_z),1);
all_rc_area_ratio_number_of_outliers = zeros(numel(unique_z),1);
all_rc_perimeter_ratio_number_of_outliers = zeros(numel(unique_z),1);
% to generate histogram counts we need to define bin edges
edges = [0.4:.02:1.7];
counts = zeros(numel(unique_z), numel(edges));
webopts = weboptions('Timeout', 60);

fprintf('Area And Perimeter Ratio Progress:');
fprintf(['\n' repmat('.',1,numel(unique_z)) '\n\n']);
parfor z_index = 1:numel(unique_z)
    % call the Renderer API to fetch tile information from target stack
    urlChar = sprintf('%s/owner/%s/project/%s/stack/%s/z/%.1f/tile-specs', ...
        rc.baseURL, rc.owner, rc.project, rc.stack,unique_z(z_index) );
    rc_data = webread(urlChar, webopts);
    
    urlChar = sprintf('%s/owner/%s/project/%s/stack/%s/z/%.1f/tile-specs', ...
        rcsource.baseURL, rcsource.owner, rcsource.project, rcsource.stack,unique_z(z_index) );
    rcsource_data = webread(urlChar, webopts);
    
    % AREA and PERIMETER process individual tiles: Calculate area and perimeter
    [is_rc_tile_in_rcsource, rcsource_tile_indices] = ismember({rc_data(:).tileId}, {rcsource_data(:).tileId});
    rc_data(rcsource_tile_indices==0)=[];
    rcsource_tile_indices = rcsource_tile_indices(is_rc_tile_in_rcsource);
    
    numel_rc_data = numel(rc_data);
    section_heights = zeros(numel_rc_data,1);
    section_widths = zeros(numel_rc_data,1);
    % check polygon area
    rc_tile_areas = zeros(numel_rc_data,1);
    rc_tile_perimeters = zeros(numel_rc_data,1);
    rc_tile_area_ratios = zeros(numel_rc_data,1);
    rc_tile_perimeter_ratios = zeros(numel_rc_data,1);
    rc_tile_ids = cell(numel_rc_data,1);
    tile_positions_transformed = cell(numel_rc_data,1);
    
    for rc_tile_index = 1:numel(rc_data)
        rc_tile = tile(rc_data(rc_tile_index));
        rcsource_tile_index = rcsource_tile_indices(rc_tile_index);
        rcsource_tile = tile(rcsource_data(rcsource_tile_index));
        rc_tile_ids{rc_tile_index} = rc_tile.renderer_id;
        % make four corners for this tile
        x = 0;
        y = 0;
        rc_tile_position_x = [x; x + rc_tile.W; x + rc_tile.W; x];
        rc_tile_position_y = [y; y    ; y + rc_tile.H; y + rc_tile.H];
        %%% transform points
        if strcmp(class(rc_tile.tform), 'affine2d')
            rc_tile_position_transformed = [rc_tile_position_x(:) rc_tile_position_y(:) [1 1 1 1]']*rc_tile.tform.T;
        else
            rc_tile_position_transformed = transformPointsInverse(rc_tile.tform,[rc_tile_position_x Py]);
        end
        tile_positions_transformed{rc_tile_index} = {rc_tile_position_transformed};
        
        section_heights(rc_tile_index) = rcsource_tile.H;
        section_widths(rc_tile_index)  = rcsource_tile.W;
        % check polygon area
        rc_tile_areas(rc_tile_index) = polyarea(rc_tile_position_transformed(:,1), rc_tile_position_transformed(:,2));
        rc_tile_area_ratios(rc_tile_index) = rc_tile_areas(rc_tile_index)/(rcsource_tile.H * rcsource_tile.W);
        %%% polygonperimeter
        rc_tile_perimeters(rc_tile_index) = rc_tile_perimeters(rc_tile_index) + sqrt((rc_tile_position_transformed(1,1)-rc_tile_position_transformed(2,1)).^2 + (rc_tile_position_transformed(1,2)-rc_tile_position_transformed(2,2)).^2);
        rc_tile_perimeters(rc_tile_index) = rc_tile_perimeters(rc_tile_index) + sqrt((rc_tile_position_transformed(2,1)-rc_tile_position_transformed(3,1)).^2 + (rc_tile_position_transformed(2,2)-rc_tile_position_transformed(3,2)).^2);
        rc_tile_perimeters(rc_tile_index) = rc_tile_perimeters(rc_tile_index) + sqrt((rc_tile_position_transformed(3,1)-rc_tile_position_transformed(4,1)).^2 + (rc_tile_position_transformed(3,2)-rc_tile_position_transformed(4,2)).^2);
        rc_tile_perimeters(rc_tile_index) = rc_tile_perimeters(rc_tile_index) + sqrt((rc_tile_position_transformed(1,1)-rc_tile_position_transformed(4,1)).^2 + (rc_tile_position_transformed(1,2)-rc_tile_position_transformed(4,2)).^2);
        rc_tile_perimeter_ratios(rc_tile_index) = rc_tile_perimeters(rc_tile_index)/(2 * rcsource_tile.H + 2* rcsource_tile.W);
    end
    all_rc_heights{z_index} = section_heights;
    all_rc_widths{z_index} = section_widths;
    counts(z_index,:) = histc(rc_tile_area_ratios, edges);
    all_rc_tile_areas{z_index} = rc_tile_areas;
    all_rc_tile_perimeters{z_index} = rc_tile_perimeters;
    all_rc_tile_area_ratios{z_index} = rc_tile_area_ratios;
    all_rc_tile_perimeter_ratios{z_index} = rc_tile_perimeter_ratios;
    [all_rc_area_ratio_median(z_index), all_rc_area_ratio_outliers_tile_ids{z_index}, all_rc_area_ratio_number_of_outliers(z_index)] = calculate_median_and_outliers( rc_tile_area_ratios, options.outlier_deviation_for_ratios, rc_tile_ids, 'fixed_cutoff');
    [all_rc_perimeter_ratio_median(z_index), all_rc_perimeter_ratio_outliers_tile_ids{z_index}, all_rc_perimeter_ratio_number_of_outliers(z_index)] = calculate_median_and_outliers( rc_tile_perimeter_ratios, options.outlier_deviation_for_ratios, rc_tile_ids, 'fixed_cutoff');
    all_rc_section_map{z_index} = tile_positions_transformed;   % needed to plot tile boxes
    fprintf('\b|\n');
end

output_struct.Area.areas = all_rc_tile_areas;
output_struct.Area.ratios = all_rc_tile_area_ratios;
output_struct.Area.outlier_tile_ids = all_rc_area_ratio_outliers_tile_ids;
output_struct.Area.number_of_outliers = all_rc_area_ratio_number_of_outliers;
output_struct.Area.median_of_means = all_rc_area_ratio_median;

output_struct.Perimeter.perimeters = all_rc_tile_perimeters;
output_struct.Perimeter.ratios = all_rc_tile_perimeter_ratios;
output_struct.Perimeter.outlier_tile_ids = all_rc_perimeter_ratio_outliers_tile_ids;
output_struct.Perimeter.number_of_outliers = all_rc_perimeter_ratio_number_of_outliers;
output_struct.Perimeter.median_of_means = all_rc_perimeter_ratio_median;

if nargout > 1, varargout{1} = all_rc_heights; end
if nargout > 2, varargout{2} = all_rc_widths; end
if nargout > 3, varargout{3} = counts; end
if nargout > 4, varargout{4} = all_rc_section_map; end
if nargout > 5, varargout{5} = section_ids_grouped_by_z; end
if nargout > 6, varargout{6} = unique_z; end

end

