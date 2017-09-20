function point_match_optimization(rc, tile_1_id, tile_2_id, SIFT_options, SURF_options, url_options, result_output_directory, figure_visibility)
%% Find point-matches for two tiles using a variety of SIFT and/or SURF parameters
% Calculates point-matches for tile_1_id and tile_2_id in rc using
% SIFT_options and SURF_options. url_options is used to render the tiles.
% An image is created with both tiles and the point matches between them.
% The figure's visibility is set to figure_visiblity and it is saved to
% result_output_directory in a subdirectory titled based on the stack and
% tile ids. Arrays of values can be used for the SIFT and SURF options
% fields to allow for testing of multiple parameter sets. eg.
% SIFT_options.renderScale = [0.4, 0.5, 0.6];
%
% SIFT_options struct fields and their defaults:
%     SIFTfdSize
%        SIFT feature descriptor size: how many samples per row and column
%        Default: 89
%     SIFTmaxScale
%        SIFT maximum scale: minSize * minScale < size < maxSize * maxScale
%        Default: 0.85
%     SIFTminScale
%        SIFT minimum scale: minSize * minScale < size < maxSize * maxScale
%        Default: 0.5
%     SIFTsteps
%        SIFT steps per scale octave
%        Default: 3
%     fillWithNoise
%        Fill each canvas image with noise before rendering to improve point match
%        derivation
%        Default: true
%     matchIterations
%        Match filter iterations
%        Default: 1000
%     matchMaxEpsilon
%        Minimal allowed transfer error for match filtering
%        Default: 20.0
%     matchMaxNumInliers
%        Maximum number of inliers for match filtering
%     matchMaxTrust
%        Reject match candidates with a cost larger than maxTrust * median cost
%        Default: 3.0
%     matchMinInlierRatio
%        Minimal ratio of inliers to candidates for match filtering
%        Default: 0.0
%     matchMinNumInliers
%        Minimal absolute number of inliers for match filtering
%        Default: 4
%     matchModelType
%        Type of model for match filtering
%        Default: AFFINE
%        Possible Values: [TRANSLATION, RIGID, SIMILARITY, AFFINE]
%     matchRod
%        Ratio of distances for matches
%        Default: 0.92
%     renderScale
%        Render canvases at this scale
%        Default: 1.0
%     outputDirectory
%        Parent directory in which subdirectories will be created to store
%        images and point-match results from SIFT.
%        NO DEFAULT VALUE
%
% SIFT_options struct fields and their defaults:
%     renderScale
%        NO DEFAULT VALUE
%
% url_options struct and their defaults:
%     normalizeForMatching
%        Default: true
%     renderWithFilter
%        Default: true
%     renderWithoutMask
%        Default: true
%     fullScaleWidth
%        Unused is default
%     fullScaleHeight         : 
%        Unused is default
% Author: David Ackerman
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Make the output directory and run SIFT and SURF if applicable
sorted_tile_ids = sort({tile_1_id, tile_2_id});
result_output_directory = [result_output_directory '/' rc.stack '/' sorted_tile_ids{1} '_and_' sorted_tile_ids{2}];
system(['mkdir -p ' result_output_directory]);
if ~isempty(SIFT_options)
    find_SIFT_point_matches(rc, tile_1_id, tile_2_id, SIFT_options, url_options, result_output_directory, figure_visibility);
end
if ~isempty(SURF_options)
    find_SURF_point_matches(rc, tile_1_id, tile_2_id, SURF_options, url_options, result_output_directory, figure_visibility);
end
end
%% Helper Functions
%% Perform all steps of SIFT from getting the tile urls to saving the final image
function find_SIFT_point_matches(rc, tile_1_id, tile_2_id, SIFT_options, url_options, result_output_directory, figure_visibility)
tile_1_url = get_tile_image_url(rc, tile_1_id, url_options, false);
tile_2_url = get_tile_image_url(rc, tile_2_id, url_options, false);

% To parallelize it, need to figure out how many total parameter sets need
% to be tested since each field may have multiple values
total_number_of_combinations = 1;
user_set_fields = fieldnames(SIFT_options);
num_fields = numel(user_set_fields);
number_of_elements_in_fields = ones(numel(user_set_fields),1);
for current_field_index = 1:num_fields
    current_field = SIFT_options.(user_set_fields{current_field_index});
    if ~(ischar(current_field) || isstring(current_field)) % if it is string or char then it is single element, otherwise it would be cell array
        total_number_of_combinations = total_number_of_combinations*numel(current_field);
        number_of_elements_in_fields(current_field_index) = numel(current_field);
    end
end
all_combinations_array = zeros(total_number_of_combinations, num_fields);
for i=1:total_number_of_combinations
    for current_field_index = 1:num_fields
        all_combinations_array(i,current_field_index) = mod(floor((i-1)/prod(number_of_elements_in_fields(current_field_index+1:end))),number_of_elements_in_fields(current_field_index))+1;
    end
end
% Loop through all the parameter set combinations
parfor i = 1:total_number_of_combinations
    SIFT_options_current = SIFT_options;
    % All of the SIFT output will be saved in a subdirectory with the name
    % of the stack
    SIFT_options_current.outputDirectory = [SIFT_options_current.outputDirectory '/' rc.stack '/' ];
    % Set SIFT_options_current based on the current parameter set
    % Also set a string of the parameter values to be used for directory
    % and image names as well as the string to be used for the legend in
    % the figure
    parameter_string = [];
    legend_string = [];
    for current_field_index = 1:num_fields
        current_field = user_set_fields{current_field_index};
        if ~strcmp(current_field, 'outputDirectory')
            parameter_string = [parameter_string (current_field) '_'];
            if iscell(SIFT_options_current.(current_field))
                SIFT_options_current.(current_field) = SIFT_options.(current_field){all_combinations_array(i,current_field_index)};
                parameter_string = [parameter_string SIFT_options_current.(current_field) '_'];
                legend_string = [legend_string, {[current_field ': ' SIFT_options_current.(current_field)]}];
            elseif ischar(SIFT_options_current.(current_field)) || isstring(SIFT_options_current.(current_field))
                SIFT_options_current.(current_field) = SIFT_options.(current_field);
                parameter_string = [parameter_string SIFT_options_current.(current_field) '_'];
                legend_string = [legend_string, {[current_field ': ' SIFT_options_current.(current_field)]}];
            else
                SIFT_options_current.(current_field) = SIFT_options.(current_field)(all_combinations_array(i,current_field_index));
                parameter_string = [parameter_string num2str(SIFT_options_current.(current_field)) '_'];
                legend_string = [legend_string, {[current_field ': ' num2str(SIFT_options_current.(current_field))]}];
            end
        end
    end
    
    % Create appropriate output directory, perform SIFT and save image if
    % point matches were found
    parameter_string(end) = [];
    SIFT_options_current.outputDirectory = [SIFT_options_current.outputDirectory parameter_string];
    system(['mkdir -p ' SIFT_options_current.outputDirectory]);
    system(['rm ' SIFT_options_current.outputDirectory '/matches.json']);
    debug_canvas_pair_matches(SIFT_options_current, tile_1_url, tile_2_url);
    fid = fopen([SIFT_options_current.outputDirectory '/matches.json']);
    im1 = rgb2gray(imread([SIFT_options_current.outputDirectory '/' tile_1_id '.jpg']));
    im2 = rgb2gray(imread([SIFT_options_current.outputDirectory '/' tile_2_id '.jpg']));
    if fid~=-1 %then it found point matches
        raw = fread(fid, inf);
        str = char(raw');
        fclose(fid);
        data = JSON.parse(str);
        pm1 = SIFT_options_current.renderScale*[cell2mat(data{1}.matches.p{1})' cell2mat(data{1}.matches.p{2})'];
        pm2 = SIFT_options_current.renderScale*[cell2mat(data{1}.matches.q{1})' cell2mat(data{1}.matches.q{2})'];
        fh = figure('visible', figure_visibility);
        showMatchedFeatures(im1, im2, pm1, pm2, 'montage');
    else % no point matches found
        pm1 = []; pm2=[];
        fh = figure('visible', figure_visibility);
        imshow(imfuse(im1,im2,'montage'));
    end
    ax = fh.CurrentAxes;
    ax.Title.String = [{'SIFT'},{['Number of Point-Matches: ' num2str(size(pm1,1))]}];
    pos_fig = get(fh, 'position');
    set(fh, 'position',  pos_fig+ [0 0 0 50]);
    annotation('textbox', [0.65, 0, 1, 0.85], 'string', legend_string, 'fitboxtotext','on','backgroundcolor','white');
    if ~isempty(result_output_directory)
        saveas(fh, [ result_output_directory '//SIFT_' parameter_string '.tif']);
    end
end
end

%% Perform SIFT on two tiles
function debug_canvas_pair_matches(SIFT_options, tile_1_url, tile_2_url)
cmd = ['/groups/flyTEM/flyTEM/render/bin/debug_canvas_pair_matches.sh ' SIFT_options.outputDirectory];
% Using these if statements allows the default values to be used when
% necessary
if isfield(SIFT_options, 'SIFTfdSize'), cmd = [cmd ' --SIFTfdSize ' num2str(SIFT_options.SIFTfdSize)]; end
if isfield(SIFT_options, 'SIFTmaxScale'), cmd = [cmd ' --SIFTmaxScale ' num2str(SIFT_options.SIFTmaxScale)]; end
if isfield(SIFT_options, 'SIFTminScale'), cmd = [cmd ' --SIFTminScale ' num2str(SIFT_options.SIFTminScale)]; end
if isfield(SIFT_options, 'SIFTsteps'), cmd = [cmd ' --SIFTsteps ' num2str(SIFT_options.SIFTsteps)]; end
if isfield(SIFT_options, 'fillWithNoise'), cmd = [cmd ' --fillWithNoise ' SIFT_options.fillWithNoise]; end
if isfield(SIFT_options, 'matchIterations'), cmd = [cmd ' --matchIterations ' num2str(SIFT_options.matchIterations)]; end
if isfield(SIFT_options, 'matchMaxEpsilon'), cmd = [cmd ' --matchMaxEpsilon ' num2str(SIFT_options.matchMaxEpsilon)]; end
if isfield(SIFT_options, 'matchMaxNumInliers'), cmd = [cmd ' --matchMaxNumInliers ' num2str(SIFT_options.matchMaxNumInliers)]; end
if isfield(SIFT_options, 'matchMaxTrust'), cmd = [cmd ' --matchMaxTrust ' num2str(SIFT_options.matchMaxTrust)]; end
if isfield(SIFT_options, 'matchMinInlierRatio'), cmd = [cmd ' --matchMinInlierRatio ' num2str(SIFT_options.matchMinInlierRatio)]; end
if isfield(SIFT_options, 'matchMinNumInliers'), cmd = [cmd ' --matchMinNumInliers ' num2str(SIFT_options.matchMinNumInliers)]; end
if isfield(SIFT_options, 'matchModelType'), cmd = [cmd ' --matchModelType ' SIFT_options.matchModelType]; end
if isfield(SIFT_options, 'matchRod'), cmd = [cmd ' --matchRod ' num2str(SIFT_options.matchRod)]; end
if isfield(SIFT_options, 'renderScale'), cmd = [cmd ' --renderScale ' num2str(SIFT_options.renderScale)]; end

cmd = [cmd ' ''' tile_1_url ''' ''' tile_2_url '''']; %adding quotes to urls
system(cmd);

end

%% Perform SURF on two tiles
function find_SURF_point_matches(rc, tile_1_id, tile_2_id, SURF_options,  url_options, result_output_directory, figure_visibility)
% get images at 100%, then only have to read them once, and scale in the
% loop
tile_1_url = get_tile_image_url(rc, tile_1_id, url_options, true);
tile_2_url = get_tile_image_url(rc, tile_2_id, url_options, true);
im1 = rgb2gray(imread(tile_1_url, 'jpg'));
im2 = rgb2gray(imread(tile_2_url, 'jpg'));
for i = 1:numel(SURF_options.renderScale)
    current_im1 = imresize(im1, SURF_options.renderScale(i));
    current_im2 = imresize(im2, SURF_options.renderScale(i));
    suppress_error = true;
    [m12_2, m12_1, fh]=find_point_matches_with_matlab(current_im1, current_im2, 'SURF', figure_visibility, suppress_error);
    parameter_string = ['renderScale_' num2str(SURF_options.renderScale(i))];
    title_string = [{['renderScale: ' num2str(SURF_options.renderScale(i))]}];
    parameter_string = strrep(parameter_string,'.','p');
    ax = fh.CurrentAxes;
    pos_fig = get(fh, 'position');
    ax.Title.String = strrep([{'SURF'}, {['Number of Point-Matches: ' num2str(size(m12_2, 1))]}], '_', '\_');
    set(fh, 'position',  pos_fig+ [0 0 0 50]);
    annotation('textbox', [0.65, 0, 1, 0.85], 'string', title_string, 'fitboxtotext','on','backgroundcolor','white');
    saveas(fh, [ result_output_directory '//SURF_' parameter_string '.tif']);
end
end