function [area_ratio_median, perimeter_ratio_median, section_map, montage_residuals_vector, tile_areas, tile_perimeters, tile_ids_vector,residuals_matrix, ...
    section_confidence, residual_outliers_tile_ids, area_ratio_outliers_tile_ids, outliers_tile_ids, Table] =...
    updated_gen_diagnostics(rcsource, rc, zstart, zend, point_matches, options)
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
    options.min_points = 5;
    options.nbrs = 4;
    options.show_deformation = 1;
    options.show_residuals = 1;
    options.show_deformation_summary = 1;
    options.show_residual_histogram = 1;
    options.nstd = 2;
    options.number_of_cross_sections;
end

%%% defaults and overrides
if ~isfield(options, 'show_residual_histogram'), options.show_residual_histogram = 0;end
if ~isfield(options, 'nstd'), options.nstd = 2;end
if ~isfield(options, 'residual_info'), options.residual_info = 0;end

% addd defaullttt forrrrrrrrrrrrrrrrrrrrrrrrrr thisssssssssssssssssssssssss
  dir_scratch = [options.dir_scratch '/temp_' num2str(randi(3000000))];
    kk_mkdir(dir_scratch);
    cd(dir_scratch);
[unique_z, section_Ids_grouped_by_z, section_Ids_all, z_all, ns1] = get_section_ids(rc, zstart, zend);

%% find corresponding tile centers
tic
% initialize variables to store deformation and residuals
section_map  = cell(numel(unique_z),1);
residual_outliers_tile_ids = cell(numel(unique_z),1);
area_ratio_outliers_tile_ids = cell(numel(unique_z),1);
outliers_tile_ids = cell(numel(unique_z),1);
section_confidence = cell(numel(unique_z), 1);
tile_areas = cell(numel(unique_z),1);
tile_perimeters = cell(numel(unique_z),1);
tile_ids_vector = cell(numel(unique_z),1);
residual_x = cell(numel(unique_z),1);
residual_y = cell(numel(unique_z),1);
tile_residuals_vector = cell(numel(unique_z),1);
montage_residuals_vector = zeros(numel(unique_z,1));
height = cell(numel(unique_z),1);
width = cell(numel(unique_z),1);

% to generate histogram counts we need to define bin edges
edges = [0.4:.02:1.7];
counts = zeros(numel(unique_z), numel(edges));
webopts = weboptions('Timeout', 60);
parfor z_index = 1:numel(unique_z)
    rc_tile_area = [];    % surface area of tiles
    tile_area_ratio = [];% surface area ratio of tiles
    tile_perimeter_ratio  = [];    % perimeter of tiles
    tile_positions_transformed = [];     % tile patch
    section_heights = [];
    section_widths = [];
    % call the Renderer API to fetch tile information from target stack
    urlChar = sprintf('%s/owner/%s/project/%s/stack/%s/z/%.1f/tile-specs', ...
        rc.baseURL, rc.owner, rc.project, rc.stack,unique_z(z_index) );
    rc_data = webread(urlChar, webopts);
    
    urlChar = sprintf('%s/owner/%s/project/%s/stack/%s/z/%.1f/tile-specs', ...
        rcsource.baseURL, rcsource.owner, rcsource.project, rcsource.stack,unique_z(z_index) );
    rcsource_data = webread(urlChar, webopts);
    
    % AREA and PERIMETER process individual tiles: Calculate area and perimeter
    [is_rc_tile_in_rcsource, rcsource_tile_indices] = ismember({rc_data(:).tileId}, {rcsource_data(:).tileId});
    rc_tiles = tile;
    rcsource_tiles = tile;
    sectionID = rc_data(1).layout.sectionId;
    count = 1;
    tile_ids = {};
% % % %             rcsource_data_for_tiles_in_rc = rcsource_data(rcsource_tile_indices);        
% % % %             x = zeros(1,numel(rc_data));
% % % %             y = zeros(1,numel(rc_data));
% % % %             rc_tile_position_x = [x; x + [rc_data.width]; x + [rc_data.width]; x];
% % % %             rc_tile_position_y = [y; y; y + [rc_data.height]; y + [rc_data.height]];
% % % %            %FOR NOW ONLY DO AFFINE%
% % % %            rc_tile_position_transformed = [rc_tile_position_x
            
    for rc_tile_index = 1:numel(rc_data)
        
        
        
        %%%%WARNING: FIND PROPER INDEX FIRST --- SOSI
        % fixed this but should I be using Msection instead?

        
        
        
        
        
        
        rc_tiles(rc_tile_index) = tile(rc_data(rc_tile_index));
        rcsource_tile_index = rcsource_tile_indices(rc_tile_index);
        rcsource_tiles(rc_tile_index) = tile(rcsource_data(rcsource_tile_index));
        tile_ids(rc_tile_index) = {rc_tiles(rc_tile_index).renderer_id};
        % make four corners for this tile
        x = 0;
        y = 0;
        rc_tile_position_x = [x; x + rc_tiles(rc_tile_index).W; x + rc_tiles(rc_tile_index).W; x];
        rc_tile_position_y = [y; y    ; y + rc_tiles(rc_tile_index).H; rc_tiles(rc_tile_index).H];
        %%% transform points
        if strcmp(class(rc_tiles(rc_tile_index).tform), 'affine2d')
            rc_tile_position_transformed = [rc_tile_position_x(:) rc_tile_position_y(:) [1 1 1 1]']*rc_tiles(rc_tile_index).tform.T;
        else
            rc_tile_position_transformed = transformPointsInverse(rc_tiles(rc_tile_index).tform,[rc_tile_position_x Py]);
        end
        tile_positions_transformed{rc_tile_index} = {rc_tile_position_transformed};
        %cm(jix,:) = sum([P(:,1)/4 P(:,2)/4],1);  % center of mass of each tile
        
        section_heights(count) = rcsource_tiles(rc_tile_index).H;
        section_widths(count)  = rcsource_tiles(rc_tile_index).W;
        % check polygon area
        rc_tile_area(count) = polyarea(rc_tile_position_transformed(:,1), rc_tile_position_transformed(:,2));
        tile_area_ratio(count) = rc_tile_area(count)/(rcsource_tiles(rc_tile_index).H * rcsource_tiles(rc_tile_index).W);
        
        %%% polygonperimeter
        rc_tile_perimeter = 0;
        rc_tile_perimeter = rc_tile_perimeter + sqrt((rc_tile_position_transformed(1,1)-rc_tile_position_transformed(2,1)).^2 + (rc_tile_position_transformed(1,2)-rc_tile_position_transformed(2,2)).^2);
        rc_tile_perimeter = rc_tile_perimeter + sqrt((rc_tile_position_transformed(2,1)-rc_tile_position_transformed(3,1)).^2 + (rc_tile_position_transformed(2,2)-rc_tile_position_transformed(3,2)).^2);
        rc_tile_perimeter = rc_tile_perimeter + sqrt((rc_tile_position_transformed(3,1)-rc_tile_position_transformed(4,1)).^2 + (rc_tile_position_transformed(3,2)-rc_tile_position_transformed(4,2)).^2);
        rc_tile_perimeter = rc_tile_perimeter + sqrt((rc_tile_position_transformed(1,1)-rc_tile_position_transformed(4,1)).^2 + (rc_tile_position_transformed(1,2)-rc_tile_position_transformed(4,2)).^2);
        tile_perimeter_ratio(count) = rc_tile_perimeter/(2 * rcsource_tiles(rc_tile_index).H + 2* rcsource_tiles(rc_tile_index).W);
        
        
        count = count + 1;
    end
    height{z_index} = section_heights;
    width{z_index} = section_widths;
    tile_ids_vector{z_index} = tile_ids;
    counts(z_index,:) = histc(tile_area_ratio, edges, 2);
    tile_areas{z_index} = rc_tile_area;
    tile_perimeters{z_index} = tile_perimeter_ratio;
    area_ratio_median(z_index) = median(tile_area_ratio);
    perimeter_ratio_median(z_index) = median(tile_perimeter_ratio);
    section_map{z_index} = tile_positions_transformed;   % needed to plot tile boxes
     %% determine area outliers
    % for each tile, calculate the mean of means of point-match residuals
    area_ratio_mean = mean(tile_area_ratio);
    area_ratio_std = std(tile_area_ratio);
        I = logical(abs(tile_area_ratio-area_ratio_mean)>options.nstd*area_ratio_std);
    %I = bsxfun(@gt, abs(bsxfun(@minus, Aratio, meann)), opts.nstd*stdd);
    area_ratio_outliers = find(I);
    area_ratio_outliers_tile_ids{z_index} = tile_ids_vector{z_index}(area_ratio_outliers);
    
    if options.residual_info
        %% %%% determine point-matches, solution and residuals for this section
        
        % First: load point-matches and section into "L" (point-matches are in L's pm struct field)
        if (z_index + options.nbrs+1)>numel(unique_z)
            point_match_zlast = unique_z(end);
        else
            point_match_zlast = unique_z(z_index + options.nbrs);
        end
        [L]  = ...
            load_point_matches(unique_z(z_index), point_match_zlast, rc, point_matches, options.nbrs, ...
            options.min_points, 0);
        
        % initialize variable transformed_point_match_residuals to hold point-match residuals
        transformed_point_match_residuals = zeros(numel(L.tiles),1);   % transformed_point_match_residuals === transformed point-match residuals
        visits = zeros(numel(L.tiles),1); % stores the number of tile-pair visits every tile experiences
        tile_residuals = cell(numel(L.tiles),1);
        % Second: generate point-match residuals from L.pm by transforming them and taking the sum of squared
        % residuals
        for point_match_index = 1:size(L.pm.M,1)
            adjacent_tile_1 = L.pm.adj(point_match_index,1);
            adjacent_tile_2 = L.pm.adj(point_match_index,2);
            point_matches_tile_1 = L.pm.M{point_match_index,1};
            point_matches_tile_2 = L.pm.M{point_match_index,2};
            point_matches_tile_1 = [point_matches_tile_1 ones(size(point_matches_tile_1,1),1)]*L.tiles(adjacent_tile_1).tform.T;  % apply transformation
            point_matches_tile_2 = [point_matches_tile_2 ones(size(point_matches_tile_2,1),1)]*L.tiles(adjacent_tile_2).tform.T;  % apply transformation
            residual = real(sum(sqrt((point_matches_tile_1(:,1)-point_matches_tile_2(:,1)).*(point_matches_tile_1(:,1)-point_matches_tile_2(:,1))  + (point_matches_tile_1(:,2)-point_matches_tile_2(:,2)).* (point_matches_tile_1(:,2)-point_matches_tile_2(:,2)))));    %%%% sum of squared residuals
            residual = residual/size(point_matches_tile_1,1);  % mean residual sum for this tile pair
            transformed_point_match_residuals(adjacent_tile_1) = transformed_point_match_residuals(adjacent_tile_1) + residual;  % add to bucket of tile adjacent_tile_1
            transformed_point_match_residuals(adjacent_tile_2) = transformed_point_match_residuals(adjacent_tile_2) + residual;  % add to bucket of tile adjacent_tile_2
            visits(adjacent_tile_1) = visits(adjacent_tile_1) + 1;              % counter for number of times tile a1 is visited
            visits(adjacent_tile_2) = visits(adjacent_tile_2) + 1;
            tile_residuals{adjacent_tile_1} = [tile_residuals{adjacent_tile_1} residual];  % aggregate residuals for tile a1
            tile_residuals{adjacent_tile_2} = [tile_residuals{adjacent_tile_2} residual];
        end
        tile_residuals_vector{z_index} = tile_residuals;  % store tile residuals for this section
        montage_residuals_vector(z_index) = nanmedian(transformed_point_match_residuals./visits);
        %% determine residual outliers
        tres = tile_residuals_vector{z_index};  % all tile residuals for section zu1(z_index)
        c = zeros(numel(tres),1);
        % for each tile, calculate the mean of means of point-match residuals
        for tix = 1:numel(tres)
            c(tix) = sum(tres{tix}(:))/numel(tres{tix}(:));
        end
        section_confidence{z_index} = c; %
        residual_mean = mean(section_confidence{z_index});
        residual_std = std(section_confidence{z_index});
        I = logical(abs(section_confidence{z_index}-residual_mean)>options.nstd*residual_std);
        %    I = bsxfun(@gt, abs(bsxfun(@minus, section_conf{z_index}, meann)), opts.nstd*stdd);
        residual_outliers = find(I);
        residual_outliers_tile_ids{z_index} = tile_ids_vector{z_index}(residual_outliers);
        %%%% outlier index
        outlier_index = unique([residual_outliers(:); area_ratio_outliers(:)]);
        outliers_tile_ids{z_index} = tile_ids_vector{z_index}(outlier_index);
    end
end
%% Cross-section residuals  
if options.residual_info
    residuals_matrix = diag(montage_residuals_vector);
    [T, map_id, tIds, z_val] = load_all_transformations(rc, unique_z, dir_scratch);
    %Convert T to cell array
    %T=[T,repmat([0,0,1],length(T),1)];
    offsets = [T(:,3),T(:,6)];
    offsets = mat2cell(offsets,ones(size(offsets,1),1),2);
    T(:,[3,6])=[];
    T=reshape(T, length(T),2,2);
    T=num2cell(T,[2,3]);
    T=cellfun(@squeeze,T,'UniformOutput',false);
    unique_z_number_of_elements = numel(unique_z);
    residuals_vector=zeros(unique_z_number_of_elements,options.number_of_cross_sections*2);
    parfor z_index=1:unique_z_number_of_elements
        new_residuals_values=zeros(1,options.number_of_cross_sections*2);
        for next_index = 1:options.number_of_cross_sections%z_index+1:z_index+2
            if z_index+next_index<=unique_z_number_of_elements
                 factor = options.xs_weight/(next_index+1);
                 disp('Loading transformations and tile/canvas ids from Renderer database.....');
                 [cross_section_point_matches, adjacency, weights, number_of_point_matches] = load_cross_section_pm(point_matches, section_Ids_grouped_by_z{z_index}, section_Ids_grouped_by_z{z_index+next_index}, ...
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
        
    end
    for z_index = 1:unique_z_number_of_elements
        for next_index = 1:options.number_of_cross_sections
            if z_index+next_index<=unique_z_number_of_elements
                residuals_matrix(z_index,z_index+next_index) = residuals_vector(z_index,next_index);
                residuals_matrix(z_index+next_index,z_index) = residuals_vector(z_index,next_index+options.number_of_cross_sections);
            end
        end
    end
end
% Resx_store = Resx;
% Resy_store = Resy;
% Resx = cell2mat(Resx);
% Resy = cell2mat(Resy);
%%
n = 0.1; % number of std for cutoff to determine outliers

%%
if options.show_deformation || options.show_residuals
    figure;plot(unique_z, montage_residuals_vector, '-o', 'LineWidth',2);title(['Point-match residuals: ' num2str(zstart) ' to ' num2str(zend)]);
    xlabel('Section number');
    ylabel('Point-match residual (sum of sum of squared distance)');
    figure;plot(unique_z, area_ratio_median, '-o', 'LineWidth',2);title(['Area ratio median per layer: ' num2str(zstart) ' to ' num2str(zend)]);
    xlabel('Section number');
    ylabel('Area ratio median (ideally one)');
    %figure;plot(zu1, mS, 'LineWidth',2);title('Perimeter median perlayer');
    
    maxres = 2.0;  % cap res at this value so that color differences can be discerned
    
    %     Resx(Resx>mean(Resx)+n*std(Resx))= mean(Resx) + n* std(Resx);
    %     Resy(Resy>mean(Resy)+n*std(Resy))= mean(Resy) + n* std(Resy);
    
    
    % figure; hist(real(Resx), 100); title([char(rc.stack) ' -- histogram for all tiles Resx']);axis tight
    % figure; hist(real(Resy), 100); title([char(rc.stack) ' -- histogram for all Resy']);axis tight
    
    %% display section tile box figures
    % at this point Resx and Resy contain residual information
    % tile_areas contains areas for individual tiles
    areas = [cell2mat(tile_areas')]'./[cell2mat(height').*cell2mat(width')]';
    areas = abs(1-areas);  % measure area devisation from 1
    %%% display histogram of tile areas for full slab
    figure;hist(areas,100);title('Deformation histogram: (ideally zero) for whole slab');axis tight
    xlim([0 1]);
    c_bounds = [min(0) max(20)];
    area_bounds = [0 2];
    colormap jet;
    
    for z_index = 1:numel(unique_z)  % loop over sections
        if options.show_deformation
            areas = tile_areas{z_index}./(height{z_index}.*width{z_index});
            % to properly see color differences we need to get rid of extreme values
            areas(areas>mean(areas)+n*std(areas))= mean(areas) + n* std(areas);
            draw_colored_boxes(section_map{z_index}, areas, area_bounds, ['Deformation : ' num2str(unique_z(z_index))]); % generate figure for y residuals
        end
        
        if options.show_residuals
            % to properly see color differences we need to get rid of extreme values
            tres = tile_residuals_vector{z_index};  % all tile residuals for section zu1(z_index)
            c = zeros(numel(tres),1);
            for tix = 1:numel(tres)
                c(tix) = sum(tres{tix}(:));
            end
            c(c>mean(c)+n*std(c))= mean(c) + n* std(c);
            draw_colored_boxes(section_map{z_index}, c, [0 max(c)], ['Residuals for: ' num2str(unique_z(z_index))]); % generate figure for y residuals
            
            %              resx = Resx_store{z_index};
            %             resx(resx>mean(resx)+n*std(resx))= mean(resx) + n* std(resx);
            %              resy = Resy_store{z_index};
            %             resy(resy>mean(resy)+n*std(resy))= mean(resy) + n* std(resy);
            %
            %resx(resx>maxres) = maxres;
            % resy(resy>maxres) = maxres;
            %resx_bounds(2) = maxres;
            %resy_bounds(2) = maxres;
            %           draw_colored_boxes(sctn_map{z_index}, resx, resx_bounds, ['Residuals x: ' num2str(zu1(z_index))]); % generate figure for y residuals
            %draw_colored_boxes(sctn_map{z_index}, resy, resy_bounds, ['Residuals y: ' num2str(zu1(z_index))]); % generate figure for y residuals
        end
    end
end

%% list section outliers and statistics per section
Table = [];
if options.residual_info
disp('Area should be close to 1.0');
MeanA = zeros(numel(unique_z),1);
MedianA = zeros(numel(unique_z),1);
MedianResidual = zeros(numel(unique_z),1);
MedianPerimeter = zeros(numel(unique_z),1);
ResidualOutliers  = zeros(numel(unique_z),1);
Outliers          = zeros(numel(unique_z), 1);
AreaOutliers  = zeros(numel(unique_z),1);

for z_index = 1:numel(unique_z)  % loop over sections
    areas = abs(tile_areas{z_index}./(height{z_index}.*width{z_index}));
    perim = abs(tile_perimeters{z_index});
    
    MeanA(z_index) = mean(areas);
    MedianA(z_index) = median(areas);
    MedianResidual(z_index) = montage_residuals_vector(z_index);
    MedianPerimeter(z_index) = median(perim);
    ResidualOutliers(z_index) = numel(residual_outliers_tile_ids{z_index});
    AreaOutliers(z_index) = numel(area_ratio_outliers_tile_ids{z_index});
    Outliers(z_index) = numel(outliers_tile_ids{z_index});
    SectionName{z_index} = num2str(unique_z(z_index));
%     str1 = [num2str(zu1(z_index)) ' Mean A: ' num2str(mean(areas)) ...
%         ' ----- Median A: ' num2str(median(areas)) ' ---- Median Residual (px): '...
%         num2str(confidence(z_index)) ' Median Perimeter: ' num2str(median(perim))];
%     disp([str1]);
    
end
Table = table(MeanA,MedianA,MedianResidual,MedianPerimeter, ResidualOutliers, AreaOutliers, Outliers,...
    'RowNames',SectionName');
disp(Table);
end
%% generate tile-based residual measure and potential outliers
% section_conf = {};
% for z_index = 1:numel(zu1)  % loop over sections
%     
%     tres = res_tiles_vec{z_index};  % all tile residuals for section zu1(z_index)
%     c = zeros(numel(tres),1);
%     % for each tile, calculate the mean of means of point-match residuals
%     for tix = 1:numel(tres)
%         c(tix) = sum(tres{tix}(:))/numel(tres{tix}(:));
%     end
%     
%     section_conf{z_index} = c;
%     
%     % filter by deviation from mean
%     %c(c>mean(c)+n*std(c))= mean(c) + n* std(c);
%     %% determine outliers
%     meann = mean(section_conf{z_index});
%     stdd = std(section_conf{z_index});
%     I = bsxfun(@gt, abs(bsxfun(@minus, section_conf{z_index}, meann)), 2*stdd);
%     outliers_tid{z_index} = tidsvec{z_index}(find(I));
%     
% endix
%%
% sc = cell2mat(section_conf);
% 
% if opts.show_residual_histogram
%     hist(sc, 100);
% end
%% summarize deformation for whole stack
if options.show_deformation_summary
    cc = counts(1:end,:);
    % plot results
    hf = figure;
    ha = axes;
    b = bar3(edges, cc.'); % note the transpose to get the colors right
    xlabel('section number')
    ylabel('Area ratio (ideally 1)');
    zlabel('Area ratio');
    
    
    for k = 1:length(b)
        zdata = b(k).ZData;
        b(k).CData = zdata;
        b(k).FaceColor = 'interp';
    end
    view(2);
    
    title('summary area deformation of stack');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% analyze results
% zconf = [];
% for z_index = 1:1%numel(zu1)
%     c = confidence{z_index};
%     t = tidsvec{z_index};
%     for tix = 1:size(c,1)
%         if isempty(c{tix})
%             zconf(tix,:) = [nan nan];
%         else
%             zconf(tix,:) = [sum(c{tix},1)/size(c{tix},1)];
%         end
%         disp([num2str(z_index) ' '  num2str(tix) ' ' t{tix} '  ' num2str(zconf(tix,:))]);
%     end
% end

%[obj, h, rh, A, minconf, maxconf] = show_map_color(obj, parm, conf,  minconf, maxconf);

%% look at section map
% c = [];
% Ao = 2560*2160;
% colormap jet;
% for z_index = 45:45 %numel(zu1)
%     %figure;
%     sm = sctn_map{z_index};
%     areas = tile_areas{z_index}/Ao;
%     %caxis([min(areas) max(areas)]);
%     c = mat2gray(areas, [min(areas) max(areas)]);
%
%     for tix = 1:numel(sm)
%         P = sm{tix}{1};   % patch for this tile
%         patch( P(:,1), P(:,2), [c(tix)],'FaceColor', 'flat',   'EdgeColor', 'k' , 'Facealpha', 0.4);
%     end
%     daspect([1 1 1]);colorbar; axis ij
% end

%%
% zu = zu1;
% % p = polyfit(zu,mA,2);
% xo = 1;
% yo = 1.05;
% n = 2;
% Aeq = xo.^(n:-1:0);
% beq = yo;
% %
% % Aeq = [];
% % beq = [];
% V = [];C = [];
% V (:,n+1) = ones(length(zu), 1, class(zu));
% for jx = n:-1:1, V(:, jx) = zu(:).*V(:, jx+1);end
% C = V;
% mA = mA(:);
% p = lsqlin(C,mA,[],[],Aeq, beq);
% fac = polyval(p, zu);
% plot(zu,mA);hold on;plot(zu, fac);
%
%
%
%
%
%
%












