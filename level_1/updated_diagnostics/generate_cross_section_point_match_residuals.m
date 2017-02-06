function [residuals_matrix] =...
    generate_cross_section_point_match_residuals(rc, zstart, zend, point_matches, options, unique_z, section_Ids_grouped_by_z)
% Calculates cross section residuals for renderer collection rc, from
% zstart to zend using point matches point_matches. Need to ensure units
% are correct.
% opts fields and their defaults:
%    min_points               : 10
%    max_points               : 100
%    number_of_cross_sections : 2
%    xs_weight                : 0.5
%    dir_scratch              : '/scratch/ackermand'     
% Output:
%       residuals_matrix
% Author: David Ackerman
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4
    % Configure point-match collection
    point_matches.server           = 'http://10.40.3.162:8080/render-ws/v1';
    point_matches.owner            = 'flyTEM';
    point_matches.match_collection = 'v12_dmesh';
    
end
if nargin<5
    % Configure options
    options.min_points = 10;
    options.max_points = 100;
    options.number_of_cross_sections=2; % This is how many sections above/below the current section we will calculate residuals between
    options.xs_weight=0.5;
    options.dir_scratch = '/scratch/ackermand';
end

if ~isfield(options, 'plot_cross_section_residuals')
   options.plot_cross_section_residuals = false;
end

if nargin<7
    % Get section ids
    [unique_z, section_Ids_grouped_by_z, ~, ~, ~] = get_section_ids(rc, zstart, zend); 
end

% Create scratch directory and change to it
 dir_scratch = [options.dir_scratch '/temp_' num2str(randi(3000000))];
 kk_mkdir(dir_scratch);
 cd(dir_scratch);

webopts = weboptions('Timeout', 60);

% Create simplified residuals matrix where the residuals for the nth
% section (Sn) will be stored as, eg.:
%  Column 1   |   Column 2   |   Column 3   |   Column 4
% Sn -> Sn+1  |  Sn -> Sn+2  |  Sn+1 -> Sn  |  Sn+2 - > Sn 
unique_z_number_of_elements = numel(unique_z);
simplified_residuals_matrix=zeros(unique_z_number_of_elements,options.number_of_cross_sections*2);

% Print out status and loop through all unique zs
fprintf('Progress:\n');
fprintf(['\n' repmat('.',1,unique_z_number_of_elements) '\n\n']);
parfor z_index=1:unique_z_number_of_elements
    % Load in transformations and map ids necessary to calculate
    % residuals for current section
    [T, map_id, ~, ~] = load_all_transformations(rc, unique_z(z_index:min(z_index+options.number_of_cross_sections,unique_z_number_of_elements)), dir_scratch);
    % new_residuals_values will store the residuals as:
    % [ Sn -> Sn+1 ,  Sn -> Sn+2  ,  Sn+1 -> Sn  ,  Sn+2 -> Sn ]
    % Necessary to have this variable because of parfor
    new_residuals_values=zeros(1,options.number_of_cross_sections*2);
    % Will want everything as cells to be able to use cellfun below.
    % We store the offsets as their own array and convert the remaining
    % part of T to its own cell array.
    offsets = [T(:,3),T(:,6)];
    offsets = mat2cell(offsets,ones(size(offsets,1),1),2);
    T(:,[3,6])=[];
    T=reshape(T, length(T),2,2);
    T=num2cell(T,[2,3]);
    T=cellfun(@squeeze,T,'UniformOutput',false);
    % Loop through all the required sections and calculate cross-sectional
    % residuals
    for section_step = 1:options.number_of_cross_sections
        if z_index+section_step<=unique_z_number_of_elements
            % Load the cross section point matches and adjacency matrix
            factor = options.xs_weight/(section_step+1);
            [cross_section_point_matches, adjacency, ~, ~] = load_cross_section_pm(point_matches, section_Ids_grouped_by_z{z_index}, section_Ids_grouped_by_z{z_index+section_step}, ...
                map_id, options.min_points, options.max_points, webopts, factor);
            if ~isempty(cross_section_point_matches)
                % We now calculate the residuals for each tile pair as the
                % average of all point match residuals between the two
                % tiles. Below, the point match positions for tile 1 are
                % pm1*t1+offsets1 and for tile2 are pm2*t2+offsets2. We
                % access the appropriate transformations and offsets using
                % the adjacency matrix.
                residuals = cellfun(@(pm1,pm2,t1,t2,offsets1,offsets2) mean( sqrt( sum(((pm1*t1+offsets1)-(pm2*t2+offsets2)).^2,2) ) ), ...
                    cross_section_point_matches(:,1), cross_section_point_matches(:,2), T(adjacency(:,1)), T(adjacency(:,2)), offsets(adjacency(:,1)), offsets(adjacency(:,2)));
            
                % To calculate the section residuals, we want to take the
                % median of the mean residuals between each tile in the
                % section and all its adjacent tiles. 
                [~,~,ic] = unique(adjacency(:,1)); % Now each repeated element of adjacency(:,1) is indexed by the same value in ic, in the range 1 to length(unique(adjacency(:,1))). If we want ic to correspond to the first time the value appears in adjacency, use 'stable'
                new_residuals_values(section_step) = median(accumarray(ic,residuals,[],@mean)); % Calculates the mean for each tile in the section (all residuals which have matching indices in ic), and calculates the median of those means
                [~,~,ic] = unique(adjacency(:,2)); % Repeat for eg Sn+1 ->  Sn, which can differ slightly from Sn -> Sn+1
                new_residuals_values(section_step+options.number_of_cross_sections) = median(accumarray(ic,residuals,[],@mean));
            end
            simplified_residuals_matrix(z_index,:) = new_residuals_values;
        end
    end
    fprintf('\b|\n');
end

% Create the full residuals matrix
residuals_matrix = zeros(length(unique_z));
for z_index = 1:unique_z_number_of_elements
    for section_step = 1:options.number_of_cross_sections
        if z_index+section_step<=unique_z_number_of_elements
            residuals_matrix(z_index, z_index+section_step) = simplified_residuals_matrix(z_index,section_step);
            residuals_matrix(z_index+section_step, z_index) = simplified_residuals_matrix(z_index,section_step+options.number_of_cross_sections);
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
    max_residuals_matrix = max(residuals_matrix(:));
    imshow(residuals_matrix/max_residuals_matrix);
    label_spacing = max(1,floor(length(all_z)/10));
    set(gca,'XTick',(1:label_spacing:length(all_z)),'XTickLabel',all_z(1:label_spacing:end));
    set(gca,'YTick',(1:label_spacing:length(all_z)),'YTickLabel',all_z(1:label_spacing:end));
    c=colorbar();
    ylabel(c,'Residuals');
    set(c,'YTick',(0:.25:1),'YTickLabel',(0:max_residuals_matrix/4:max_residuals_matrix));

end

end











