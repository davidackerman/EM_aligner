% % % diagnostics
clear all;
options.dir_scratch = '/scratch/ackermand';
options.xs_weight = 0.5;
options.min_points            = 8;
options.max_points            = Inf;
options.number_of_cross_sections = 0;
options.show_residuals = false; % Show_residuals and show_deformations are for plotting section maps
options.show_deformations = false;
options.show_table = false; % Whether or not to show the table
options.output_data_per_tile = true; % Store all the data per tile (eg. residuals, area, area ratio etc) or not
options.outlier_deviation_for_residuals = 10; % Cutoff average residual for tile, beyond which it is considered to be an outlier
options.outlier_deviation_for_ratios = 0.10; % Cutoff for area ratio and perimeter ratio outliers: any tile ratios that stray by more than outlier_deviation_for_ratios*100% are outliers
options.visible = false;
options.save_comparison_figure = true;
options.save_comparison_figure_directory = '/groups/flyTEM/home/ackermand/beautification/temp_david2/images/';
options.save_comparison_text_directory = '/groups/flyTEM/home/ackermand/beautification/temp_david2/text/';
options.verbose = false;

% Source, used for area and perimeter ratios, empty if don't want area and
% perimeter calculations
rcsource.baseURL = 'http://10.37.5.60:8080/render-ws/v1';
rcsource.owner = 'flyTEM';
rcsource.project = 'FAFB00';
rcsource.stack = 'v12_acquire_merged';
sl.source_collection = rcsource;

rcs(1).owner = 'flyTEM';
rcs(1).project = 'FAFB00';%'FAFB00';
rcs(1).stack = 'v12_acquire_merged';%'v13_montage';
rcs(1).baseURL = 'http://10.37.5.60:8080/render-ws/v1';
rcs(1).verbose = 0;

rcs(2).owner = 'flyTEM';
rcs(2).project = 'FAFBv14_kk';%'FAFB00';
rcs(2).stack = 'kk_13_pm7';%'v13_montage';
rcs(2).baseURL = 'http://10.37.5.60:8080/render-ws/v1';
rcs(2).verbose = 0;

rcs(3).owner = 'flyTEM';
rcs(3).project = 'FAFBv14_kk';%'FAFB00';
rcs(3).stack = 'kk14_montage';%'v13_montage';
rcs(3).baseURL = 'http://10.37.5.60:8080/render-ws/v1';
rcs(3).verbose = 0;



% rcs(2).owner = 'flyTEM';
% rcs(2).project = 'FAFB00_beautification';%'FAFB00';
% rcs(2).stack = 'Revised_FAFB_montage_remove_small_clusters';%'v13_montage';
% rcs(2).baseURL = 'http://10.37.5.60:8080/render-ws/v1';
% rcs(2).verbose = 0;

% Point matches
pm.server = 'http://10.40.3.162:8080/render-ws/v1';
pm.owner            = 'flyTEM';
pm.match_collection = 'v12_dmesh';
pm.verbose = 0;

 options.zstart=100;
 options.zend=105;
%%
sl.options = options;
sl.renderer_collections = rcs;
sl.point_match_collection = pm;
sl.verbose = true;

fn = [pwd '/sample_input_diagnostics_maps_and_residuals_function.json'];
str = savejson('', sl);
fid = fopen(fn,'w');
fprintf(fid,str);
fclose(fid);
%% make sure we can read this file
loaded = loadjson(fileread(fn));
diagnostics_maps_and_residuals_function(fn);
