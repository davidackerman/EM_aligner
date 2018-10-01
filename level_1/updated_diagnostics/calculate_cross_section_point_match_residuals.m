function [median_residuals_matrix, max_residuals_matrix] = calculate_cross_section_point_match_residuals(rc, point_matches, zstart, zend, varargin)
%% Generate statistics about cross section residuals
% Calculates cross section residuals for renderer collection rc, using
% zstart, zend, point matches point matches and varargin. If varargin does
% not contain unique_z and section_ids_grouped_by_z, they will be obtained
% based on zstart and zend. 
% varargin should contain 0-3 arguments:
%    1:                       The input should be options
%    2:                       The input should be unique_z and section_ids_grouped_by_z
%    3:                       The input should be unique_z, section_ids_grouped_by_z and options
% opts fields and their defaults:
%    min_points                   : 10
%    max_points                   : 100
%    number_of_cross_sections     : 2
%    xs_weight                    : 0.5
%    dir_scratch                  : '/scratch/ackermand' 
%    plot_cross_section_residuals : false
%    filter_point_matches         : true
%    verbose                      : true
% Output:
%    residuals_matrix
% Author: David Ackerman
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check the input arguments
options = [];
if isempty(varargin)
    [unique_z, section_ids_grouped_by_z, ~, ~, ~] = get_section_ids(rc, zstart, zend+1);
    z_too_large = (unique_z>=zend+1);
    unique_z(z_too_large) = [];
    section_ids_grouped_by_z(z_too_large) = [];
elseif length(varargin)==1
    if isstruct(varargin{1})
        options = varargin{1};
        [unique_z, section_ids_grouped_by_z, ~, ~, ~] = get_section_ids(rc, zstart, zend+1);
        z_too_large = (unique_z>=zend+1);
        unique_z(z_too_large) = [];
        section_ids_grouped_by_z(z_too_large) = [];
    else
        error('Not enough input arguments: need both unique_z and section_ids_grouped_by_z.');
    end
elseif length(varargin)==2
    if isstruct(varargin{2})
        error('Not enough input arguments: need both unique_z and section_ids_grouped_by_z.');
    else
       unique_z = varargin{1};
       section_ids_grouped_by_z = varargin{2};
    end
elseif length(varargin)==3
    if isstruct(varargin{3})
        unique_z = varargin{1};
        section_ids_grouped_by_z = varargin{2};
        options = varargin{3};
    end
end

% Group zs when necessary
floor_unique_z = floor(unique_z);
unique_merged_z = unique(floor_unique_z);
section_ids_grouped_by_z_and_section_id = section_ids_grouped_by_z;
for i=1:numel(section_ids_grouped_by_z_and_section_id)
    for j=i+1:numel(section_ids_grouped_by_z_and_section_id)
        if floor(unique_z(i)) == floor(unique_z(j)) %Then two zs share the same "z"
            unique_sections = ~ismember(section_ids_grouped_by_z_and_section_id{j}, section_ids_grouped_by_z_and_section_id{i});
            section_ids_grouped_by_z_and_section_id{i}(end+1:end+sum(unique_sections)) = section_ids_grouped_by_z_and_section_id{j}(unique_sections);
            section_ids_grouped_by_z_and_section_id{j}(:) = [];
        end
    end
end
section_ids_grouped_by_z_and_section_id = section_ids_grouped_by_z_and_section_id(~cellfun('isempty',section_ids_grouped_by_z_and_section_id));

% Configure options
new_dir_scratch = false;
if ~isfield(options, 'min_points'), options.min_points = 10; end
if ~isfield(options, 'max_points'), options.max_points = 100; end
if ~isfield(options, 'number_of_cross_sections'), options.number_of_cross_sections = 2; end % This is how many sections above/below the current section we will calculate residuals between
if ~isfield(options, 'xs_weight'), options.xs_weight = 0.5; end
if ~isfield(options, 'dir_scratch')
    new_dir_scratch=true;
    options.dir_scratch = [pwd '/scratch_' num2str(randi(10000)) '_' datestr(datetime('now'),'yyyymmdd_HHMMSS')];
    warning('Will create temporary scratch directory %s which will be cleaned after', options.dir_scratch);
end
if ~isfield(options, 'plot_cross_section_residuals'), options.plot_cross_section_residuals = false; end
if ~isfield(options, 'filter_point_matches'), options.filter_point_matches = true; end
if ~isfield(options, 'verbose'), options.verbose = true; end
if ~isfield(options, 'pmopts')
    options.pmopts.NumRandomSamplingsMethod = 'Desired confidence';
    options.pmopts.MaximumRandomSamples = 3000;
    options.pmopts.DesiredConfidence = 99.5;
    options.pmopts.PixelDistanceThreshold = 1;
end
if ~isfield(options, 'do_cross_only') || ~options.do_cross_only, options.do_cross_only = true; end 
if ~isfield(options, 'nbrs') || options.nbrs~=1, options.nbrs = 1; end % only want to consider one other section at a time

% Create scratch directory and change to it
 dir_current = pwd;
 dir_scratch = [options.dir_scratch '/temp_' num2str(randi(3000000))];
 kk_mkdir(dir_scratch);
 cd(dir_scratch);

webopts = weboptions('Timeout', 60);

% Create simplified residuals matrix where the residuals for the nth
% section (Sn) will be stored as, eg.:
%  Column 1   |   Column 2   |   Column 3   |   Column 4
% Sn -> Sn+1  |  Sn -> Sn+2  |  Sn+1 -> Sn  |  Sn+2 - > Sn 
unique_merged_z_number_of_elements = numel(unique_merged_z);
simplified_median_residuals_matrix=zeros(unique_merged_z_number_of_elements,options.number_of_cross_sections*2);
simplified_max_residuals_matrix=zeros(unique_merged_z_number_of_elements,options.number_of_cross_sections*2);

% Print out status and loop through all unique zs
if options.verbose
    fprintf('Cross Section Residuals Progress:');
    fprintf(['\n' repmat('.',1,unique_merged_z_number_of_elements) '\n\n']);
end
for z_index=1:unique_merged_z_number_of_elements
    % Load in transformations and map ids necessary to calculate
    % residuals for current section
    cross_section_z_max = unique_merged_z(z_index)+options.number_of_cross_sections;
    
    [T, map_id, tIds, z_val, r, c] = load_all_transformations(rc, unique_z(floor_unique_z>= unique_merged_z(z_index) & floor_unique_z<=cross_section_z_max), dir_scratch);
    
    % new_residuals_values will store the residuals as:
    % [ Sn -> Sn+1 ,  Sn -> Sn+2  ,  Sn+1 -> Sn  ,  Sn+2 -> Sn ]
    % Necessary to have this variable because of parfor
    new_median_residuals_values=NaN(1,options.number_of_cross_sections*2);
    new_max_residuals_values=NaN(1,options.number_of_cross_sections*2);
    % Store the offsets of T separately, and convert the remainder to a
    % cell array for efficiency
    offsets = [T(:,3),T(:,6)];
    T(:,[3,6])=[];
    T=reshape(T, size(T,1),2,2);
    T=num2cell(T,[2,3]);
    T=cellfun(@squeeze,T,'UniformOutput',false);
    % Loop through all the required sections and calculate cross-sectional
    % residuals
    for section_step = 1:options.number_of_cross_sections
        if z_index+section_step<=unique_merged_z_number_of_elements && (unique_merged_z(z_index+section_step) - unique_merged_z(z_index) == section_step)
            % Load the cross section point matches and adjacency matrix and filter           
            factor = options.xs_weight/(section_step+1);
            zs = [unique_merged_z(z_index), unique_merged_z(z_index+section_step)];
            sIDs = [section_ids_grouped_by_z_and_section_id(z_index), section_ids_grouped_by_z_and_section_id(z_index+section_step)];
            PM=[];
            [PM.M, PM.adj, PM.W, PM.np] = system_solve_helper_load_point_matches(zs, options, point_matches, map_id, sIDs, size(T,1), r, c);
            PMold = [];
            [PMold.M, PMold.adj, PMold.W, PMold.np] = load_cross_section_pm(point_matches, section_ids_grouped_by_z_and_section_id{z_index}, section_ids_grouped_by_z_and_section_id{z_index+section_step}, ...
                map_id, options.min_points, options.max_points, webopts, factor);
            if options.filter_point_matches
                if isfield(options, 'pmopts')
                    PMold = filter_pm(PMold, options.pmopts);
                else
                    PMold = filter_pm(PMold);
                end
            end
            %PM.M is cross section point match residuals, PM.adj is
            %adjacency matrix
            if ~isempty(PM.M)
                % We now calculate the residuals for each tile pair as the
                % average of all point match residuals between the two
                % tiles. Below, we store the information for each pair for,
                % eg., Sn->Sn+1 in ..._current_section_tiles and Sn+1->Sn
                % in ..._next_section_tiles. ic1 assigns each unique tile
                % number in Sn an index between 1 and
                % length(unique(adjacency(:,1)). Likewise for ic2.
                [~,~,ic1] = unique(PM.adj(:,1),'stable');
                num_unique_current_section_tiles = max(ic1);
                [~,~,ic2] = unique(PM.adj(:,2),'stable');
                num_unique_next_section_tiles = max(ic2);
                residuals_current_section_tiles = zeros(num_unique_current_section_tiles,1);
                counts_current_section_tiles = zeros(num_unique_current_section_tiles,1);
                residuals_next_section_tiles = zeros(num_unique_next_section_tiles,1);
                counts_next_section_tiles = zeros(num_unique_next_section_tiles,1);
                for i=1:size(PM.M,1)
                    residual = mean(sqrt(sum((PM.M{i,1}*T{PM.adj(i,1)}+offsets(PM.adj(i,1),:) - (PM.M{i,2}*T{PM.adj(i,2)}+offsets(PM.adj(i,2),:))).^2,2)));
                    residuals_current_section_tiles(ic1(i)) = residuals_current_section_tiles(ic1(i)) + residual;
                    residuals_next_section_tiles(ic2(i)) = residuals_next_section_tiles(ic2(i)) + residual;
                    counts_current_section_tiles(ic1(i)) = counts_current_section_tiles(ic1(i))+1;
                    counts_next_section_tiles(ic2(i)) = counts_next_section_tiles(ic2(i))+1;
                end
                
                % To calculate the section residuals, we want to take the
                % median of the mean residuals; residuals/counts provides
                % the mean residuals between each tile and all its adjacent
                % tiles, as stored above
                new_median_residuals_values(section_step) = median(residuals_current_section_tiles./counts_current_section_tiles);
                new_median_residuals_values(section_step+options.number_of_cross_sections) = median(residuals_next_section_tiles./counts_next_section_tiles);
                new_max_residuals_values(section_step) = max(residuals_current_section_tiles./counts_current_section_tiles);
                new_max_residuals_values(section_step+options.number_of_cross_sections) = max(residuals_next_section_tiles./counts_next_section_tiles);
            end
            simplified_median_residuals_matrix(z_index,:) = new_median_residuals_values;
            simplified_max_residuals_matrix(z_index,:) = new_max_residuals_values;
        end
    end
    if options.verbose, fprintf('\b|\n'); end
end

% Create the full residuals matrix
median_residuals_matrix = NaN(length(unique_merged_z));
max_residuals_matrix = NaN(length(unique_merged_z));
for z_index = 1:unique_merged_z_number_of_elements
    for section_step = 1:options.number_of_cross_sections
        if z_index+section_step<=unique_merged_z_number_of_elements
            median_residuals_matrix(z_index, z_index+section_step) = simplified_median_residuals_matrix(z_index,section_step);
            median_residuals_matrix(z_index+section_step, z_index) = simplified_median_residuals_matrix(z_index,section_step+options.number_of_cross_sections);
            
            max_residuals_matrix(z_index, z_index+section_step) = simplified_max_residuals_matrix(z_index,section_step);
            max_residuals_matrix(z_index+section_step, z_index) = simplified_max_residuals_matrix(z_index,section_step+options.number_of_cross_sections);
        end
    end
end

% Plot residuals
if options.plot_cross_section_residuals
    all_z = (zstart:zend);
    iptsetpref('ImshowAxesVisible','on')
    figure();
    Outer = get(gca,'OuterPosition');
    set(gca,'OuterPosition',[Outer(1)+0.05 Outer(2) + 0.05 Outer(3)-0.05 Outer(4)-0.05])
    max_residuals_matrix = max(median_residuals_matrix(:));
    imshow(median_residuals_matrix/max_residuals_matrix);
    label_spacing = max(1,floor(length(all_z)/10));
    set(gca,'XTick',(1:label_spacing:length(all_z)),'XTickLabel',all_z(1:label_spacing:end));
    set(gca,'YTick',(1:label_spacing:length(all_z)),'YTickLabel',all_z(1:label_spacing:end));
    c=colorbar();
    ylabel(c,'Residuals');
    set(c,'YTick',(0:.25:1),'YTickLabel',(0:max_residuals_matrix/4:max_residuals_matrix));
    title(['Cross Section Residuals For ' num2str(zstart) '-' num2str(zend)]); 
end
cd(dir_current);
if new_dir_scratch
    system(sprintf('rm -rf %s', options.dir_scratch));
end
end











