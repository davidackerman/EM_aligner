function create_section_diagnostic_maps( rc_in, output_struct_in, section_zs, options, section_maps_in )
if numel(rc_in) ~= numel(output_struct_in), error('rc and output_struct need to be the same size'); end
if nargin<5, section_maps_provided = false; else section_maps_provided = true; end
if section_maps_provided
    if numel(rc_in) ~= numel(section_maps_in)
        error('rc and section_maps need to be the same size');
    end
end
if ~isfield(options, 'save_comparison_figure'), options.save_comparison_figure = false; end
if ~isfield(options, 'store_residual_outlier_pair_information'), options.store_residual_outlier_pair_information = false; end

num_to_compare = numel(rc_in);
if isfield(output_struct_in(1), 'Area')
    number_of_rows = 3;
    number_of_columns = num_to_compare;
else
    number_of_rows = num_to_compare;
    number_of_columns = 1;
end
all_residuals=cell(1,numel(rc_in));
all_outlier_count=zeros(1,numel(rc_in));
for z_index=1:numel(section_zs)
    if isfield(options, 'visible') && options.visible == true
        new_fig = figure();
    else
        new_fig = figure('visible','off');
    end
    for index_to_compare = 1:numel(rc_in)
        rc = rc_in(index_to_compare);
        output_struct = output_struct_in(index_to_compare);
        map_options=[];
        current_z = section_zs(z_index);
        if section_maps_provided
            all_section_maps = section_maps_in(index_to_compare);
        else
            all_section_maps = calculate_section_maps(rc, current_z, current_z);
        end
        map_options.subplotting = true;
        plot_count = index_to_compare;
        subplot(number_of_rows,number_of_columns,plot_count);
        
        if isfield(output_struct, 'Area')
            ratio_bounds = [1-options.outlier_deviation_for_ratios,1+options.outlier_deviation_for_ratios];
            areas = output_struct.Area.ratios{z_index};
            areas(areas<=ratio_bounds(1)) = -Inf;
            areas(areas>=ratio_bounds(2)) = Inf;
            outlier_count_small = sum(areas<=ratio_bounds(1));
            outlier_count_large = sum(areas>=ratio_bounds(2));
            map_options.label_str{1} = [{['\fontsize{12}' strrep(rc.stack,'_','\_')]} ...
                {['\fontsize{11}Area Ratios for z: ' num2str(current_z)] }...
                {sprintf('\\color[rgb]{0 0 0} %d Outliers \\leq %.2f \\color[rgb]{1 0 0} %d Outliers \\geq %.2f',outlier_count_small, ratio_bounds(1), outlier_count_large, ratio_bounds(2)) }];
            map_options.label_str{2} = 'Area Ratio';
            draw_colored_boxes(rc, current_z, all_section_maps{z_index}, areas, ratio_bounds, map_options);
            
            plot_count = plot_count+number_of_columns;
            subplot(number_of_rows,number_of_columns,plot_count);
            perimeters = output_struct.Perimeter.ratios{z_index};
            perimeters(perimeters<=ratio_bounds(1)) = -Inf;
            perimeters(perimeters>=ratio_bounds(2)) = Inf;
            outlier_count_small = sum(perimeters<=ratio_bounds(1));
            outlier_count_large = sum(perimeters>=ratio_bounds(2));
            map_options.label_str{1} = [{['Perimeter Ratios for z: ' num2str(current_z)] }...
                {sprintf('\\color[rgb]{0 0 0} %d Outliers \\leq %.2f \\color[rgb]{1 0 0} %d Outliers \\geq %.2f',outlier_count_small, ratio_bounds(1), outlier_count_large, ratio_bounds(2)) }];
            map_options.label_str{2} = 'Perimeter Ratio';
            draw_colored_boxes(rc, current_z, all_section_maps{z_index}, perimeters, ratio_bounds, map_options);
            plot_count = plot_count+number_of_columns;
        end
        subplot(number_of_rows,number_of_columns,plot_count);
        max_residuals = cellfun(@max,output_struct.MontageResiduals.values{z_index},'UniformOutput',false);  % all tile residuals for section zu1(z_index)
        empties = cellfun('isempty',max_residuals);
        max_residuals(empties) = {NaN};
        max_residuals = cell2mat(max_residuals);
        residuals=max_residuals(output_struct.MontageResiduals.all_mapping_from_section_map_indices_to_current_indicies{z_index});
        residuals_bounds = [0, options.outlier_deviation_for_residuals];%[min(residuals), options.outlier_deviation_for_residuals];
        outlier_count_large = sum(residuals>=residuals_bounds(2));
        map_options.only_greater_than = true;
        map_options.label_str{1} = [{['Residuals for z: ' num2str(current_z)]}...
            {sprintf('\\color[rgb]{0 0 0} %d Unconnected Tiles \\color[rgb]{1 0 0} %d Outliers \\geq %0.2f ',output_struct.MontageResiduals.unconnected_count(z_index), outlier_count_large,residuals_bounds(2))} ];
        if ~isfield(output_struct, 'Area')
            map_options.label_str{1} = [{['\fontsize{12}' strrep(rc.stack,'_','\_')]} map_options.label_str{1}];
        end
        map_options.label_str{2} = 'Residuals (Pixels)';
        draw_colored_boxes(rc, current_z, all_section_maps{z_index}, residuals, residuals_bounds, map_options); % generate figure for y residuals
        all_residuals{index_to_compare}=residuals;
        all_outlier_count(index_to_compare) = outlier_count_large;
    end
end
pos=get(new_fig, 'position');
pos(3) = 922;
pos(4) = 1031;
set(new_fig, 'position', pos);
if options.save_comparison_figure
    if ~exist(options.save_comparison_text_directory), mkdir(options.save_comparison_text_directory); end
    if ~isfield(options, 'save_comparison_figure_directory'), error('Need save_comparison_figure_directory if saving'); end
    if ~strcmp(options.save_comparison_figure_directory(end),'/'), options.save_comparison_figure_directory(end+1) = '/'; end
    if ~exist(options.save_comparison_figure_directory), mkdir(options.save_comparison_figure_directory); end
    %if numel(rc_in) > 1
    output_base_string='';
    for rc_index=1:numel(rc_in)
        output_base_string=[output_base_string rc_in(rc_index).stack '_vs_'];
    end
    output_base_string=output_base_string(1:end-4);
    
    saveas(gcf,[options.save_comparison_figure_directory '/' output_base_string '_section_' num2str(section_zs(z_index)) '.tif']);
    fileid=fopen([options.save_comparison_text_directory '/' output_base_string '_section_' num2str(section_zs(z_index)) '.txt'],'w');
    for rc_index=1:numel(rc_in)
        fprintf(fileid,'%d %d %d %f %f %f %f %f \n', numel(all_residuals{rc_index}), numnan(all_residuals{rc_index}), all_outlier_count(rc_index), numnan(all_residuals{rc_index})*1.0/numel(all_residuals{rc_index}), all_outlier_count(rc_index)*1.0/numel(all_residuals{rc_index}), nanmean(all_residuals{rc_index}), nanmedian(all_residuals{rc_index}), nanvar(all_residuals{rc_index}));
    end
    fclose(fileid);
    
    if options.store_residual_outlier_pair_information
        for index_to_compare = 1:numel(rc_in)
            rc = rc_in(index_to_compare);
            outlier_pm_max_pair_data = output_struct_in(index_to_compare).MontageResiduals.outlier_pm_max_pair_data;
            
            for z_index=1:numel(outlier_pm_max_pair_data)
                if ~isempty(outlier_pm_max_pair_data{z_index})
                    fileid = fopen([options.save_comparison_text_directory '/' rc.stack '_outlier_pm_pair_data_section_' num2str(section_zs(z_index)) '.txt'],'w');
                    fprintf(fileid, '%s \n', outlier_pm_max_pair_data{z_index}{:});
                    fclose(fileid);
                end
            end
        end
    end
        
%     else
%         saveas(gcf,[options.save_comparison_figure_directory '/' rc_in(1).stack '_section_' num2str(section_zs(z_index)) '.tif']);
%         fileid=fopen([options.save_comparison_text_directory '/' rc_in(1).stack '_section_' num2str(section_zs(z_index)) '.txt'],'w');
%         fprintf(fileid,'%d %d %d %f %f %f %f %f \n', numel(all_residuals{:}), numnan(all_residuals{:}), all_outlier_count, numnan(all_residuals{:})*1.0/numel(all_residuals{:}), all_outlier_count*1.0/numel(all_residuals{:}), nanmean(all_residuals{:}), nanmedian(all_residuals{:}), nanvar(all_residuals{:}));
%         fclose(fileid);
%     end
end
if ~options.visible, close(new_fig); end
end


