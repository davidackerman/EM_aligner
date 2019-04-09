function diagnostics_maps_and_residuals_function(varargin)
if nargin==1
    inputs = loadjson(fileread(varargin{1}));
    options=inputs.options;
    rcs_temp=inputs.renderer_collections;
    if numel(rcs_temp)>1
        for i=1:numel(rcs_temp)
            rcs(i)=rcs_temp{i};
        end
    else
        rcs=rcs_temp;
    end
        
    pm=inputs.point_match_collection;
    if isfield(inputs,'source_collection')
        rcsource=inputs.source_collection;
    else
        rcsource=[];
    end
    if inputs.verbose,
        disp('Running diagonstics');
        disp(['Using input file: ' varargin{1}]);
        disp('Using options:');disp(options);
        disp('Using collections:')
        for i=1:numel(rcs)
            disp(rcs(i));
        end
        if isfield(inputs,'source_collection')
            disp('Using source collection:');disp(rcsource);
        end
    end
else
    options=varargin{1};
    rcs=varargin{2}; 
    pm=varargin{3};
end

if nargin==4
    rcsource=varargin{4};
elseif nargin==3
    rcsource=[];
end

if isfield(options,'zstart')
    [zu, section_ids_grouped_by_z, ~, ~, ~] = get_section_ids(rcs(1), options.zstart, options.zend);
elseif isfield(options,'zs')
    zu=[];
    section_ids_grouped_by_z=[];
    for z_current=options.zs
        [zu_current, section_ids_grouped_by_z_current, ~, ~, ~] = get_section_ids(rcs(1), z_current, z_current);
        zu=[zu, zu_current];
        section_ids_grouped_by_z=[section_ids_grouped_by_z, section_ids_grouped_by_z_current];
    end
end
if ~isfield(options, 'store_residual_outlier_pair_information'), options.store_residual_outlier_pair_information = false; end
%% beautified
%options.save_comparison_figure_directory = '/groups/flyTEM/home/ackermand/beautification/temp';%v14_montages/diagnostics_presentation';
numel_zu=numel(zu);
    
parfor z_index = 1:numel_zu
    z = zu(z_index);
    [output_structs, section_maps] = updated_gen_diagnostics(rcsource, rcs(1), pm, z, z, options);
    for rc_index=2:numel(rcs)
        [output_structs(rc_index), section_maps(rc_index)] = updated_gen_diagnostics(rcsource, rcs(rc_index), pm, z, z, options);
    end
    create_section_diagnostic_maps(rcs, output_structs, z,options, section_maps)
end

sections=zeros(numel_zu,1);
is_merged=zeros(numel_zu,1);
all_num_tiles=zeros(numel_zu,numel(rcs));
all_num_unconnected_tiles=zeros(numel_zu,numel(rcs));
all_num_residual_outliers=zeros(numel_zu,numel(rcs));
all_fraction_unconnected_tiles=zeros(numel_zu,numel(rcs));
all_fraction_residual_outliers=zeros(numel_zu,numel(rcs));
all_mean_residuals=zeros(numel_zu,numel(rcs));
all_median_residuals=zeros(numel_zu,numel(rcs));
all_variance_residuals=zeros(numel_zu,numel(rcs));
output_base_string='';
for rc_index=1:numel(rcs)
    output_base_string=[output_base_string rcs(rc_index).stack '_vs_'];
end
output_base_string=output_base_string(1:end-4);
parfor z_index = 1:numel(zu)
    z = zu(z_index);
    file = dir([options.save_comparison_text_directory '/' output_base_string '_section_' num2str(z) '.txt']);
    [all_num_tiles(z_index,:), all_num_unconnected_tiles(z_index,:), all_num_residual_outliers(z_index,:), all_fraction_unconnected_tiles(z_index,:), all_fraction_residual_outliers(z_index,:), all_mean_residuals(z_index,:), all_median_residuals(z_index,:), all_variance_residuals(z_index,:)] = textread([options.save_comparison_text_directory '/' file.name], '%d %d %d %f %f %f %f %f');
    sections(z_index)=z;
    is_merged(z_index)=numel(section_ids_grouped_by_z{zu==z})>1;
end

output_labels = {'Z,', 'IsMerged,'};
output_data = [sections, is_merged];
file_name='';
for i=1:numel(rcs)
    file_name=[file_name, rcs(i).stack '_vs_'];
    output_labels={output_labels{:}, [rcs(i).stack ': # Tiles,'], [rcs(i).stack ': # Unconnected,'], [rcs(i).stack ': # Residual Outliers,'], [rcs(i).stack ': Fraction Unconnected,'], [rcs(i).stack ': Fraction Residual Outliers,'], [rcs(i).stack ': Residual Mean,'], [rcs(i).stack ': Residual Median,'], [rcs(i).stack ': Residual Variance,']};
    output_data = [output_data, all_num_tiles(:,i), all_num_unconnected_tiles(:,i), all_num_residual_outliers(:,i), all_fraction_unconnected_tiles(:,i), all_fraction_residual_outliers(:,i), all_mean_residuals(:,i), all_median_residuals(:,i), all_variance_residuals(:,i)];
end
file_name(end-3:end)=[];
output_labels=cell2mat(output_labels);
output_labels(end)=[]; %remove last comma
delete([options.save_comparison_text_directory '/' file_name '_summary.csv']);
fid = fopen([options.save_comparison_text_directory '/' file_name '_summary.csv'],'w');
fprintf(fid,'%s\n',output_labels);
fclose(fid);
dlmwrite([options.save_comparison_text_directory '/' file_name '_summary.csv'],output_data,'precision',6,'-append');

if options.store_residual_outlier_pair_information
    for rc_index=1:numel(rcs) %output file for each rc
        summary_file_name = [options.save_comparison_text_directory '/' rcs(rc_index).stack '_outlier_pm_pair_data.csv'];
        delete(summary_file_name);
        summary_fid = fopen(summary_file_name,'w');
        fprintf(summary_fid,'%s\n',cell2mat({'Z,', 'Max PM Residual,', 'Tile 1,', 'Tile 2,', 'PM X,', 'PM Y,', 'Mean Tile Pair Residual Is Outlier'}));
        fclose(summary_fid);
        fid = fopen(summary_file_name,'a');
        for z = zu
            file_name = [options.save_comparison_text_directory '/' rcs(rc_index).stack '_outlier_pm_pair_data_section_' num2str(z) '.txt'];
            if exist(file_name, 'file') == 2 %otherwise there were no outliers for that section
                [string_z, string_residual, string_tile_1, string_tile_2, string_pm_x, string_pm_y, string_mean_tile_pair_residual_is_outlier] = textread(file_name, '%s %s %s %s %s %s %s');
                output_data = [string_z, string_residual, string_tile_1, string_tile_2, string_pm_x, string_pm_y, string_mean_tile_pair_residual_is_outlier];
                [nrows, ~]=size(output_data);
                for row=1:nrows
                    fprintf(fid, '%s,%s,%s,%s,%s,%s,%s\n', output_data{row,:});
                end
                delete(file_name);
            end
        end
        fclose(fid);
    end
end