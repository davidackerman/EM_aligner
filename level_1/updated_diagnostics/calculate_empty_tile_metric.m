function fft2_results  = calculate_empty_tile_metric( rc, zu, opts )
%% Determines if a tile is blank based on its radially averaged fft2 and a cutoff.
% opts has fields displayed below as: fieldname (default)
% scale (0.25):                 scale of image used for fft2
% angle_to_skip (15):           angle within x or y axis to skip during radial averaging
% min_radius (100):             minimum radius for radial averaging
% max_radius (200):             maximum radius for radial averaging
% cutoff (1E4):                 cutoff for blank tile metric for the tile to be considered blank
% save_removed_tile_images (0): whether to save images, though this doesn't work well in parfor
% output_directory (pwd):       output directory for saving images
% Returns fft2_results

if ~isfield(opts, 'scale') opts.scale = 0.25; end
if ~isfield(opts, 'angle_to_skip') opts.angle_to_skip = 15; end
if ~isfield(opts, 'min_radius') opts.min_radius = 100; end
if ~isfield(opts, 'max_radius') opts.max_radius = 200; end
if ~isfield(opts, 'cutoff_value'), opts.cutoff_value = 1E4; end
if ~isfield(opts, 'save_removed_tile_images'), opts.save_removed_tile_images = 0;  end
if ~isfield(opts, 'output_directory'), opts.output_directory = pwd; end
if ~isfield(opts, 'debug'), opts.debug = 0; end  
if isempty(zu)
   zu = get_section_ids(rc); 
end
fft2_results = calculate_all_fft2_and_profile(rc, zu, opts);%calculate_fft2_and_profile(rc, tile_ids, 0.25);
%%
if opts.debug
    all_tile_ids = [fft2_results(:).tile_ids];
    all_tile_metrics = [fft2_results(:).blank_metric];
    tile_indices = 1:numel(all_tile_metrics);
    is_in_middle = [fft2_results(:).is_in_middle];
    [sorted_metrics,sorted_indices] = sort(all_tile_metrics);
    sorted_ids = all_tile_ids(sorted_indices);
    sorted_is_in_middle = is_in_middle(sorted_indices);
    figure(); hold on;
    plot(tile_indices(sorted_is_in_middle), sorted_metrics(sorted_is_in_middle),'ro')
    plot(tile_indices(~sorted_is_in_middle), sorted_metrics(~sorted_is_in_middle),'b.')
    sorted_ids_in_middle = sorted_ids(sorted_is_in_middle);
    sorted_metrics_in_middle = sorted_metrics(sorted_is_in_middle);
    xlabel('Sorted Tile Number');
    ylabel('Blank Tile Metric');
    legend('Middle Tiles','Edge Tiles');
end
end
%%
function fft2_results =calculate_all_fft2_and_profile(rc, zs, opts)
[~, ~, tile_ids, section_ids, ~, ~, tile_image_paths] = load_all_transformations(rc(1), zs, opts.dir_scratch);
num_tiles = numel(tile_ids);
scale = opts.scale;
if(rc.verbose), fprintf(['\n' repmat('.',1,numel(zs)) '\n\n']); end
empty_metric = zeros(size(tile_ids));
if opts.save_removed_tile_images
    system(['mkdir -p ' opts.output_directory]);
end
tic;
parfor i=1:num_tiles
    img = imresize(imread(tile_image_paths{i}),scale);
    [~, empty_metric(i)] = radial_average_of_fft2(img, opts);
    if opts.save_removed_tile_images
        if empty_metric(i) < opts.cutoff_value
            prefix = sprintf('%011.4f', round(empty_metric(i)*10000)/10000);
            my_parfor_imwrite(imresize(img,0.15), [output_directory prefix '_' tile_ids{i} '.tif']);
        end
    end
    % Get average
    if(rc.verbose), fprintf('\b|\n'); end
end
toc;
fft2_results.tile_ids = tile_ids;
fft2_results.empty_metric = empty_metric;
fft2_results.z = section_ids;
end
%%
function [radial_profile, blank_metric] = radial_average_of_fft2(img, opts)
angle_to_skip = opts.angle_to_skip;
min_radius = opts.min_radius;
max_radius = opts.max_radius;
fft2_magnitude = abs(fftshift(fft2(img)));
[rows, columns] = size(fft2_magnitude);
midRow = rows/2+1;
midCol = columns/2+1;
radial_profile = zeros((max_radius-min_radius)+1, 1);
count = zeros((max_radius-min_radius)+1, 1);
for col = floor(midCol-max_radius):ceil(midCol+max_radius)
    for row = floor(midRow-max_radius):ceil(midRow+max_radius)
        radius = sqrt((row - midRow) ^ 2 + (col - midCol) ^ 2);
        if radius>=min_radius && radius<=max_radius
            current_angle_from_center = atand((row-midRow)/(col-midCol));
            if current_angle_from_center > angle_to_skip && current_angle_from_center < 90-angle_to_skip
                thisIndex = ceil(radius)-min_radius+1 ;
                radial_profile(thisIndex) = radial_profile(thisIndex) + fft2_magnitude(row, col);
                count(thisIndex) = count(thisIndex) + 1;
            end
        end
    end
end
%figure();
%imshow(log(fft2_magnitude)/10.0)
%imagesc(log(fft2_magnitude));
% Get average
radial_profile = radial_profile ./ count;
blank_metric = nanmean(radial_profile(10:30));
end

function my_parfor_imwrite(img, file_name)
imwrite(img, file_name);

end

