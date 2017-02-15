function [output_struct] = updated_gen_diagnostics(rcsource, rc, zstart, zend, point_matches, options)
%% Generates information about tile deformation and residuals
% Summarizes tile deformation and point-match residuals per tile and
% section using rcsource, r, zstart, zend, point matches and options.
% Output is a struct that contains residuals, area ratios and perimeter
% ratios, as well as their outliers. The struct also contains the residuals
% matrix for all sections.
% options fields and their defaults:
%    outlier_deviation_for_ratios     : 0.04
%    nstd            : 2
%    residual_info   : false
%    show_deformation: 1      0 = don'e show
%                             1 = display visible figure
%                             2 = save image of invisible figure to disk (not implemented yet)
%                             3 = captures invisible figure to memory (not implemented yet)
%    show_residuals: 1        0 = don't show
%                             1 = displays a visible figure
%                             2 = save image of invisiblle figure to disk (not implemented yet)
%                             3 = captures invisible figure to memory (not implemented yet)
% Output:
%       output_struct
% Author: Khaled Khairy, David Ackerman
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<6
    options=[];
end

if ~isfield(options, 'outlier_deviation_for_ratios'), options.outlier_deviation_for_ratios = 0.04; end
if ~isfield(options, 'nstd'), options.nstd = 2; end
if ~isfield(options, 'residual_info'), options.residual_info = false; end
if ~isfield(options,'show_deformation'), options.show_deformation = 0; end
if ~isfield(options,'show_residuals'), options.show_residuals = 0; end

[ output_struct, counts, all_section_map, unique_z, section_ids_grouped_by_z] = calculate_area_and_perimeter_ratios(rcsource, rc, zstart, zend, options);
%% Montage and Cross-section residuals
[ output_struct2 ] = calculate_montage_point_match_residuals(rc, point_matches,unique_z, options);
output_struct.MontageResiduals = output_struct2;
if zstart ~= zend
    residuals_matrix = calculate_cross_section_point_match_residuals(rc, zstart, zend, point_matches, unique_z, section_ids_grouped_by_z, options);
    residuals_matrix = residuals_matrix + diag(output_struct.MontageResiduals.median_of_means);
    output_struct.CrossSectionAndMontageResidualsMatrix = residuals_matrix;
end
%%
if options.show_deformation || options.show_residuals
    figure;plot(unique_z, output_struct.MontageResiduals.median_of_means, '-o', 'LineWidth',2);title(['Point-match residuals: ' num2str(zstart) ' to ' num2str(zend)]);
    xlabel('Section number');
    ylabel('Point-match residual (median of mean distance per layer)');
    figure;plot(unique_z, output_struct.Area.median, '-o', 'LineWidth',2);title(['Area ratio median per layer: ' num2str(zstart) ' to ' num2str(zend)]);
    xlabel('Section number');
    ylabel('Area ratio median (ideally one)');
    
    %% display section tile box figures
    areas = output_struct.Area.median;
    areas = abs(1-areas);  % measure area devisation from 1
    %%% display histogram of tile areas for full slab
    figure;hist(areas,100);title('Deformation histogram: (ideally zero) for whole slab');axis tight
    xlim([0 1]);
    area_bounds = [1-options.outlier_deviation_for_ratios,1+options.outlier_deviation_for_ratios];
    for z_index = 1:numel(unique_z)  % loop over sections
        if options.show_deformation
            areas = output_struct.Area.ratios{z_index};
            areas(areas<=area_bounds(1)) = -Inf;
            areas(areas>=area_bounds(2)) = Inf;
            outlier_count_small = sum(areas<=area_bounds(1));
            outlier_count_large = sum(areas>=area_bounds(2));
            draw_colored_boxes(rc, unique_z(z_index), all_section_map{z_index}, areas, area_bounds, [{['Deformation : ' num2str(unique_z(z_index))] }...
                {sprintf('\\color[rgb]{0 0 0} %d Outliers \\leq %.2f \\color[rgb]{1 0 0} %d Outliers \\geq %.2f',outlier_count_small, area_bounds(1), outlier_count_large, area_bounds(2)) }]);
        end
        
        if options.show_residuals
            residuals = cellfun(@mean,output_struct.MontageResiduals.values{z_index});  % all tile residuals for section zu1(z_index)
            residuals_mean = nanmean(residuals); 
            residuals_std = nanstd(residuals);
            residuals_bounds = [min(residuals), residuals_mean+options.nstd*residuals_std];
            outlier_count_large = sum(residuals>=residuals_bounds(2));
            only_greater_than = true;                      
            draw_colored_boxes(rc, unique_z(z_index), all_section_map{z_index}, residuals, residuals_bounds, [{['Residuals for: ' num2str(unique_z(z_index))]}...
                {sprintf('\\color[rgb]{0 0 0} %d Unconnected Tiles \\color[rgb]{1 0 0} %d Outliers \\geq %0.2f ',output_struct.MontageResiduals.unconnected_count, outlier_count_large,residuals_bounds(2))} ],... 
                only_greater_than); % generate figure for y residuals
        end
    end
end

%% list section outliers and statistics per section
Table = [];
if options.residual_info
    MedianAreaRatio = output_struct.Area.median;
    MedianPerimeterRatio = output_struct.Perimeter.median;
    MedianMontageResidual = output_struct.MontageResiduals.median_of_means;
    AreaOutliers  = output_struct.Area.outlier_count;
    PerimeterOutliers  = output_struct.Perimeter.outlier_count;
    MontageResidualOutliers  = output_struct.MontageResiduals.outlier_count;
    TotalAreaAndPerimeterOutliers = cellfun(@(area_outlier_tile_ids, perimeter_outlier_tile_ids) numel(unique([area_outlier_tile_ids(:);perimeter_outlier_tile_ids(:)])), output_struct.Area.outlier_tile_ids,  output_struct.Perimeter.outlier_tile_ids);
    SectionName = cellstr(num2str(unique_z'));
    Table = table(MedianAreaRatio,MedianPerimeterRatio,MedianMontageResidual, AreaOutliers, PerimeterOutliers, MontageResidualOutliers, TotalAreaAndPerimeterOutliers,...
        'RowNames',SectionName');
    output_struct.Table = Table;
    disp(Table);
end

%% summarize deformation for whole stack
edges = [0:.02:10];
if options.show_deformation
    if numel(unique_z) == 1
        [non_zero_rows,~] = find(counts'~=0);
        bar(edges+0.01,counts);
        xlim([edges(min(non_zero_rows)) edges(max(non_zero_rows))]);
        xlabel('Area ratio (ideally 1)')
        ylabel('Count');
        title(['Area ratio for: ' num2str(unique_z)]);
    else
        [non_zero_rows,~] = find(counts'~=0);
        cc = counts(1:end,:);
        % plot results
        hf = figure;
        ha = axes;
        b = bar3(edges+0.01, cc.', 1); % note the transpose to get the colors right
        xlabel('section number');
        ylabel('Area ratio (ideally 1)');
        zlabel('Area ratio');
        ylim([edges(min(non_zero_rows)) edges(max(non_zero_rows))]);
        colorbar
        for k = 1:length(b)
            set(b,'edgecolor','none');
            zdata = b(k).ZData;
            b(k).CData = zdata;
            b(k).FaceColor = 'interp';
        end
        view(2);
         xticklabels = get(gca,'xticklabels');
         for i=1:length(xticklabels)
            xticklabels{i} = str2num( num2str(xticklabels{i} + unique_z(1))); 
         end
         set(gca,'xticklabels', xticklabels);
         pbaspect([1 2 1])
        title('summary area deformation of stack');
    end
end











