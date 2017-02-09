% % % diagnostics
clear all;
dopts.nbrs = 0;
dopts.min_points = 3;
dopts.show_deformation_summary = 0;
dopts.dir_scratch = '/scratch/ackermand';
dopts.xs_weight = 0.5;
dopts.min_points            = 10;
dopts.max_points            = 100;
dopts.number_of_cross_sections = 2;
dopts.show_residuals = 0;
dopts.show_deformation = 0;
dopts.residual_info = 1;
dopts.nstd = 2;
dopts.outlier_deviation_for_ratios = 0.04; % Cutoff for area ratio and perimeter ratio outliers: any tile ratios that stray by more than outlier_deviation_for_ratios*100% are outliers

% Source, used for area and perimeter ratios
rcsource.baseURL = 'http://10.37.5.60:8080/render-ws/v1';
rcsource.owner = 'flyTEM';
rcsource.project = 'FAFB00';
rcsource.stack = 'v12_acquire_merged';

% Original, unbeautified stack
rc_original.baseURL = 'http://10.37.5.60:8080/render-ws/v1';
rc_original.owner = 'flyTEM';
rc_original.project = 'FAFB00';
rc_original.stack = 'v13_align';
rc_original.verbose = 0;

% First and last sections to be analyzed
nfirst = 4480; nlast = 4490;

% Beautified stack
rc_beautified = rc_original;
rc_beautified.project = 'FAFB00_beautification';
rc_beautified.stack = 'Revised_slab_4480_4490_fine';

% Point matches
pm(1).server = 'http://10.40.3.162:8080/render-ws/v1';
pm(1).owner            = 'flyTEM';
pm(1).match_collection = 'FAFB_pm_2';
pm(1).verbose = false;

pm(2).server           = 'http://10.40.3.162:8080/render-ws/v1';
pm(2).owner            = 'flyTEM';
pm(2).match_collection = 'Beautification_cross_sift_00';

original_output_struct = updated_gen_diagnostics(rcsource, rc_original, nfirst, nlast, pm, dopts);
beautified_output_struct = updated_gen_diagnostics(rcsource, rc_beautified, nfirst, nlast, pm, dopts);

%% Plot the matrices and ratio of the two

all_z = (nfirst:nlast);
iptsetpref('ImshowAxesVisible','on')
maximum_of_fine_vs_original = max([beautified_output_struct.CrossSectionAndMontageResidualsMatrix(:); original_output_struct.CrossSectionAndMontageResidualsMatrix(:)]);
label_spacing = max(1,floor(length(all_z)/5));
figure();
set(gca,'position',[0 0 1 1],'units','normalized')
subplot(3,1,1);
imshow(original_output_struct.CrossSectionAndMontageResidualsMatrix/maximum_of_fine_vs_original);
title('Original');
set(gca,'XTick',(1:label_spacing:length(all_z)),'XTickLabel',all_z(1:label_spacing:end));
set(gca,'YTick',(1:label_spacing:length(all_z)),'YTickLabel',all_z(1:label_spacing:end));
c=colorbar();
ylabel(c,'Residuals');
set(c,'YTick',(0:.25:1),'YTickLabel',(0:maximum_of_fine_vs_original/4:maximum_of_fine_vs_original));

subplot(3,1,2);
imshow(beautified_output_struct.CrossSectionAndMontageResidualsMatrix/maximum_of_fine_vs_original);
title('Beautified');
set(gca,'XTick',(1:label_spacing:length(all_z)),'XTickLabel',all_z(1:label_spacing:end));
set(gca,'YTick',(1:label_spacing:length(all_z)),'YTickLabel',all_z(1:label_spacing:end));
c=colorbar();
ylabel(c,'Residuals')
set(c,'YTick',(0:.25:1),'YTickLabel',(0:maximum_of_fine_vs_original/4:maximum_of_fine_vs_original));

subplot(3,1,3);
residuals_ratio = beautified_output_struct.CrossSectionAndMontageResidualsMatrix./original_output_struct.CrossSectionAndMontageResidualsMatrix;
max_residuals_ratio = max(residuals_ratio(:));
min_residuals_ratio = min(residuals_ratio(:));
residuals_ratio(isnan(residuals_ratio))=1;
residuals_ratio = (residuals_ratio-min_residuals_ratio)/(max_residuals_ratio-min_residuals_ratio);
index_of_unity_ratio = (1-min_residuals_ratio)*255/(max_residuals_ratio-min_residuals_ratio)+1;
red_and_blue_color = ((1:index_of_unity_ratio)-1)/(index_of_unity_ratio-1);
increased_residuals_color = ones(length(red_and_blue_color),3);
increased_residuals_color(:,[1,3]) = [red_and_blue_color', red_and_blue_color'];
red_color = flip((256:-1:index_of_unity_ratio)-index_of_unity_ratio)/(256-index_of_unity_ratio);
decreased_residuals_color = zeros(length(red_color),3);
decreased_residuals_color(:,1) = red_color; 
combined_color_map = [increased_residuals_color; decreased_residuals_color];
% combined_color_map = zeros(256,3);
% combined_color_map(end-(length(red_color)-1):end,1) = red_color;
% combined_color_map(1:length(green_color),2) = green_color;
imshow(residuals_ratio, 'colormap',combined_color_map);
title('Ratio (Beautified/Original)');
set(gca,'XTick',(1:label_spacing:length(all_z)),'XTickLabel',all_z(1:label_spacing:end));
set(gca,'YTick',(1:label_spacing:length(all_z)),'YTickLabel',all_z(1:label_spacing:end));
c=colorbar();
ylabel(c,'Ratio')
set(c,'YTick',sort([(0:.25:1),(1-min_residuals_ratio)/(max_residuals_ratio-min_residuals_ratio)]),'YTickLabel',sort([(min_residuals_ratio:(max_residuals_ratio-min_residuals_ratio)/4:max_residuals_ratio),1]));

