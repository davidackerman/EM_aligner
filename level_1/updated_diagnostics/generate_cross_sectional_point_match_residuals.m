function [residuals_matrix] =...
    generate_cross_sectional_point_match_residuals(rc, zstart, zend, point_matches, options)
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
if nargin<4
    % configure point-match collection
    point_matches.server           = 'http://10.40.3.162:8080/render-ws/v1';
    point_matches.owner            = 'flyTEM';
    point_matches.match_collection = 'v12_dmesh';
    
end
if nargin<5
    options.min_points = 10;
    options.max_points = 100;
    options.number_of_cross_sections=2;
    options.xs_weight=0.5;
    options.dir_scratch = '/scratch/ackermand';
end

%%% defaults and overrides
% if ~isfield(options, 'show_residual_histogram'), options.show_residual_histogram = 0;end
% if ~isfield(options, 'nstd'), options.nstd = 2;end
% if ~isfield(options, 'residual_info'), options.residual_info = 0;end

% addd defaullttt forrrrrrrrrrrrrrrrrrrrrrrrrr thisssssssssssssssssssssssss
dir_scratch = [options.dir_scratch '/temp_' num2str(randi(3000000))];
kk_mkdir(dir_scratch);
cd(dir_scratch);
tic;
[unique_z, section_Ids_grouped_by_z, ~, ~, ~] = get_section_ids(rc, zstart, zend);
toc;
webopts = weboptions('Timeout', 60);
%% Cross-section residuals
tic;
% if options.residual_info
%     if method_type==1
%         residuals_matrix = zeros(length(unique_z));
%         tic;
%         [T, map_id, ~, ~] = load_all_transformations(rc, unique_z, dir_scratch);
%         toc;
%         %Convert T to cell array
%         %T=[T,repmat([0,0,1],length(T),1)];
%         offsets = [T(:,3),T(:,6)];
%         offsets = mat2cell(offsets,ones(size(offsets,1),1),2);
%         T(:,[3,6])=[];
%         T=reshape(T, length(T),2,2);
%         T=num2cell(T,[2,3]);
%         T=cellfun(@squeeze,T,'UniformOutput',false);
%         unique_z_number_of_elements = numel(unique_z);
%         residuals_vector=zeros(unique_z_number_of_elements,options.options.number_of_cross_sections*2);
%         parfor z_index=1:unique_z_number_of_elements
%             new_residuals_values=zeros(1,options.options.number_of_cross_sections*2);
%             for next_index = 1:options.options.number_of_cross_sections%z_index+1:z_index+2
%                 if z_index+next_index<=unique_z_number_of_elements
%                     factor = options.xs_weight/(next_index+1);
%                     %      disp('Loading transformations and tile/canvas ids from Renderer database.....');
%                     [cross_section_point_matches, adjacency, ~, ~] = load_cross_section_pm(point_matches, section_Ids_grouped_by_z{z_index}, section_Ids_grouped_by_z{z_index+next_index}, ...
%                         map_id, options.min_points, options.max_points, webopts, factor);
%                     if ~isempty(cross_section_point_matches)
%                         residuals = cellfun(@(pm1,pm2,t1,t2,offsets1,offsets2) mean( sqrt( sum(((pm1*t1+offsets1)-(pm2*t2+offsets2)).^2,2) ) ), ...
%                             cross_section_point_matches(:,1), cross_section_point_matches(:,2), T(adjacency(:,1)), T(adjacency(:,2)), offsets(adjacency(:,1)), offsets(adjacency(:,2)));
%                         [~,~,ic] = unique(adjacency(:,1));
%                         new_residuals_values(next_index) = median(accumarray(ic,residuals,[],@mean));
%                         [~,~,ic] = unique(adjacency(:,2));
%                         new_residuals_values(next_index+options.options.number_of_cross_sections) = median(accumarray(ic,residuals,[],@mean));
%                     end
%                     residuals_vector(z_index,:) = new_residuals_values;
%                 end
%             end
%             if mod(z_index,round(unique_z_number_of_elements/100.0))==0
%                 disp('%f Percent Done\n',z_index/unique_z_number_of_elements*100);
%             end
%
%         end
%         for z_index = 1:unique_z_number_of_elements
%             for next_index = 1:options.options.number_of_cross_sections
%                 if z_index+next_index<=unique_z_number_of_elements
%                     residuals_matrix(z_index,z_index+next_index) = residuals_vector(z_index,next_index);
%                     residuals_matrix(z_index+next_index,z_index) = residuals_vector(z_index,next_index+options.options.number_of_cross_sections);
%                 end
%             end
%         end
%
%     else
residuals_matrix = zeros(length(unique_z));
unique_z_number_of_elements = numel(unique_z);
residuals_vector=zeros(unique_z_number_of_elements,options.number_of_cross_sections*2);
fprintf('Progress:\n');
fprintf(['\n' repmat('.',1,unique_z_number_of_elements) '\n\n']);
parfor z_index=1:unique_z_number_of_elements
    [T, map_id, ~, ~] = load_all_transformations(rc, unique_z(z_index:min(z_index+2,unique_z_number_of_elements)), dir_scratch);
    new_residuals_values=zeros(1,options.number_of_cross_sections*2);
    offsets = [T(:,3),T(:,6)];
    offsets = mat2cell(offsets,ones(size(offsets,1),1),2);
    T(:,[3,6])=[];
    T=reshape(T, length(T),2,2);
    T=num2cell(T,[2,3]);
    T=cellfun(@squeeze,T,'UniformOutput',false);
    for next_index = 1:options.number_of_cross_sections%z_index+1:z_index+2
        if z_index+next_index<=unique_z_number_of_elements
            factor = options.xs_weight/(next_index+1);
            %      disp('Loading transformations and tile/canvas ids from Renderer database.....');
            [cross_section_point_matches, adjacency, ~, ~] = load_cross_section_pm(point_matches, section_Ids_grouped_by_z{z_index}, section_Ids_grouped_by_z{z_index+next_index}, ...
                map_id, options.min_points, options.max_points, webopts, factor);
            if ~isempty(cross_section_point_matches)
                residuals = cellfun(@(pm1,pm2,t1,t2,offsets1,offsets2) mean( sqrt( sum(((pm1*t1+offsets1)-(pm2*t2+offsets2)).^2,2) ) ), ...
                    cross_section_point_matches(:,1), cross_section_point_matches(:,2), T(adjacency(:,1)), T(adjacency(:,2)), offsets(adjacency(:,1)), offsets(adjacency(:,2)));
                [~,~,ic] = unique(adjacency(:,1));
                new_residuals_values(next_index) = median(accumarray(ic,residuals,[],@mean));
                [~,~,ic] = unique(adjacency(:,2));
                new_residuals_values(next_index+options.number_of_cross_sections) = median(accumarray(ic,residuals,[],@mean));
            end
            residuals_vector(z_index,:) = new_residuals_values;
        end
    end
    fprintf('\b|\n');
    %   if mod(count_done,round(unique_z_number_of_elements/10.0))==0
    %       fprintf('%f Percent Done\n',count_done*100.0/unique_z_number_of_elements);
    %   end
    
end
for z_index = 1:unique_z_number_of_elements
    for next_index = 1:options.number_of_cross_sections
        if z_index+next_index<=unique_z_number_of_elements
            residuals_matrix(z_index,z_index+next_index) = residuals_vector(z_index,next_index);
            residuals_matrix(z_index+next_index,z_index) = residuals_vector(z_index,next_index+options.number_of_cross_sections);
        end
    end
end
toc;
end

%end













