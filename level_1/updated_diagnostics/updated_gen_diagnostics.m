function [area_ratio_median, perimeter_ratio_median, section_map, confidence, tile_areas, tile_perimeters, tile_ids_vector,...
    section_confidence, residual_outliers_tile_ids, area_ratio_outliers_tile_ids, outliers_tile_ids,  Table] =...
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
% Author: Khaled Khairy
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
end

%%% defaults and overrides
if ~isfield(options, 'show_residual_histogram'), options.show_residual_histogram = 0;end
if ~isfield(options, 'nstd'), options.nstd = 2;end

[unique_z, sID1, sectionId1, z1, ns1] = get_section_ids(rc, zstart, zend);

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
confidence = zeros(numel(unique_z,1));
height = cell(numel(unique_z),1);
width = cell(numel(unique_z),1);

% to generate histogram counts we need to define bin edges
edges = [0.4:.02:1.7];
counts = zeros(numel(unique_z), numel(edges));
webopts = weboptions('Timeout', 60);
for zix = 1:numel(unique_z)
    rc_tile_area = [];    % surface area of tiles
    tile_area_ratio = [];% surface area ratio of tiles
    tile_perimeter_ratio  = [];    % perimeter of tiles
    tile_positions_transformed = [];     % tile patch
    height = [];
    width = [];
    % call the Renderer API to fetch tile information from target stack
    urlChar = sprintf('%s/owner/%s/project/%s/stack/%s/z/%.1f/tile-specs', ...
        rc.baseURL, rc.owner, rc.project, rc.stack,unique_z(zix) );
    rc_data = webread(urlChar, webopts);
    
    urlChar = sprintf('%s/owner/%s/project/%s/stack/%s/z/%.1f/tile-specs', ...
        rcsource.baseURL, rcsource.owner, rcsource.project, rcsource.stack,unique_z(zix) );
    rcsource_data = webread(urlChar, webopts);
    
    % AREA and PERIMETER process individual tiles: Calculate area and perimeter
    rc_tiles = tile;
    rcsource_tiles = tile;
    sectionID = rc_data(1).layout.sectionId;
    count = 1;
    tile_ids = {};
    for rc_tile_index = 1:numel(rc_data)
        rc_tiles(rc_tile_index) = tile(rc_data(rc_tile_index));
        rcsource_tiles(rc_tile_index) = tile(rcsource_data(rc_tile_index));
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
        
        height(count) = rcsource_tiles(rc_tile_index).H;
        width(count)  = rcsource_tiles(rc_tile_index).W;
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
    height{zix} = height;
    width{zix} = width;
    tile_ids_vector{zix} = tile_ids;
    counts(zix,:) = histc(tile_area_ratio, edges, 2);
    tile_areas{zix} = rc_tile_area;
    tile_perimeters{zix} = tile_perimeter_ratio;
    area_ratio_median(zix) = median(tile_area_ratio);
    perimeter_ratio_median(zix) = median(tile_perimeter_ratio);
    section_map{zix} = tile_positions_transformed;   % needed to plot tile boxes
     %% determine area outliers
    % for each tile, calculate the mean of means of point-match residuals
    area_ratio_mean = mean(tile_area_ratio);
    area_ratio_std = std(tile_area_ratio);
        I = logical(abs(tile_area_ratio-area_ratio_mean)>options.nstd*area_ratio_std);
    %I = bsxfun(@gt, abs(bsxfun(@minus, Aratio, meann)), opts.nstd*stdd);
    area_ratio_outliers = find(I);
    area_ratio_outliers_tile_ids{zix} = tile_ids_vector{zix}(area_ratio_outliers);
    
    %% %%% determine point-matches, solution and residuals for this section
    
    % First: load point-matches and section into "L" (point-matches are in L's pm struct field)
    %options.nbrs is number of neighbors to check
    if (zix + options.nbrs+1)>numel(unique_z)
        point_match_zlast = unique_z(end);
    else
        point_match_zlast = unique_z(zix + options.nbrs);
    end
    [L]  = ...
        load_point_matches(unique_z(zix), point_match_zlast, rc, point_matches, options.nbrs, ...
        options.min_points, 0);
    
    % initialize variable transformed_point_match_residuals to hold point-match residuals
    transformed_point_match_residuals = zeros(numel(L.tiles),1);   % transformed_point_match_residuals === transformed point-match residuals
    visits = zeros(numel(L.tiles),1); % stores the number of tile-pair visits every tile experiences
    tile_residuals = cell(numel(L.tiles),1);
    % Second: generate point-match residuals from L.pm by transforming them and taking the sum of squared
    % residuals
    residual_vector_x = [];
    residual_vector_y = [];
    for point_match_index = 1:size(L.pm.M,1)
        adjacent_tile_1 = L.pm.adj(point_match_index,1);
        adjacent_tile_2 = L.pm.adj(point_match_index,2);
        point_matches_tile_1 = L.pm.M{point_match_index,1};
        point_matches_tile_2 = L.pm.M{point_match_index,2};
        point_matches_tile_1 = [point_matches_tile_1 ones(size(point_matches_tile_1,1),1)]*L.tiles(adjacent_tile_1).tform.T;  % apply transformation
        point_matches_tile_2 = [point_matches_tile_2 ones(size(point_matches_tile_2,1),1)]*L.tiles(adjacent_tile_2).tform.T;  % apply transformation
        residual = real(sqrt(sum((point_matches_tile_1(:,1)-point_matches_tile_2(:,1)).*(point_matches_tile_1(:,1)-point_matches_tile_2(:,1))  + (point_matches_tile_1(:,2)-point_matches_tile_2(:,2)).* (point_matches_tile_1(:,2)-point_matches_tile_2(:,2)))));    %%%% sum of squared residuals
        residual = residual/size(point_matches_tile_1,1);  % mean residual sum for this tile pair
        %disp(res);
        %res_vec(pmix,:) = sqrt((m1(1)-m2(1))*(m1(1)-m2(1)) + (m1(1)-m2(1))*(m1(2)-m2(2)));
        
        %%%% sosi
        %         res_vecx(pmix,:) = abs(m1(1)-m2(1));
        %         res_vecy(pmix,:) = abs(m1(2)-m2(2));
        %
        transformed_point_match_residuals(adjacent_tile_1) = [transformed_point_match_residuals(adjacent_tile_1) + residual];  % add to bucket of tile adjacent_tile_1
        transformed_point_match_residuals(adjacent_tile_2) = [transformed_point_match_residuals(adjacent_tile_2) + residual];  % add to bucket of tile adjacent_tile_2
        visits(adjacent_tile_1) = visits(adjacent_tile_1) + 1;              % counter for number of times tile a1 is visited
        visits(adjacent_tile_2) = visits(adjacent_tile_2) + 1;
        tile_residuals{adjacent_tile_1} = [tile_residuals{adjacent_tile_1} residual];  % aggregate residuals for tile a1
        tile_residuals{adjacent_tile_2} = [tile_residuals{adjacent_tile_2} residual];
    end
    tile_residuals_vector{zix} = tile_residuals;  % store tile residuals for this section
    residual_x{zix} = residual_vector_x;
    residual_y{zix} = residual_vector_y;
    % store errors/confidence in cell array "confidence"
    %confidence(zix) = 0;
    %for tix = 1:numel(L.tiles)
    confidence(zix) = median(transformed_point_match_residuals./visits);%/visits(tix);% median of residuals of tiles in this section
    %end
    %confidence(zix) = confidence(zix)/numel(L.tiles); %mean residual for tiles
     %% determine residual outliers
     tres = tile_residuals_vector{zix};  % all tile residuals for section zu1(zix)
    c = zeros(numel(tres),1);
    % for each tile, calculate the mean of means of point-match residuals
    for tix = 1:numel(tres)
        c(tix) = sum(tres{tix}(:))/numel(tres{tix}(:));
    end
    section_confidence{zix} = c; %
    area_ratio_mean = mean(section_confidence{zix});
    area_ratio_std = std(section_confidence{zix});
    I = logical(abs(section_confidence{zix}-area_ratio_mean)>options.nstd*area_ratio_std);
%    I = bsxfun(@gt, abs(bsxfun(@minus, section_conf{zix}, meann)), opts.nstd*stdd);
    resid_out = find(I);
    residual_outliers_tile_ids{zix} = tile_ids_vector{zix}(resid_out);
    %%%% outlier index
    outlier_index = unique([resid_out(:); area_ratio_outliers(:)]);
    outliers_tile_ids{zix} = tile_ids_vector{zix}(outlier_index);
end
% Resx_store = Resx;
% Resy_store = Resy;
% Resx = cell2mat(Resx);
% Resy = cell2mat(Resy);
%%
n = 0.1; % number of std for cutoff to determine outliers

%%
if options.show_deformation || options.show_residuals
    figure;plot(unique_z, confidence, '-o', 'LineWidth',2);title(['Point-match residuals: ' num2str(zstart) ' to ' num2str(zend)]);
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
    areas = abs(1-areas);  % measure area deviation from 1
    %%% display histogram of tile areas for full slab
    figure;hist(areas,100);title('Deformation histogram: (ideally zero) for whole slab');axis tight
    xlim([0 1]);
    c_bounds = [min(0) max(20)];
    area_bounds = [0 2];
    colormap jet;
    
    for zix = 1:numel(unique_z)  % loop over sections
        if options.show_deformation
            areas = tile_areas{zix}./(height{zix}.*width{zix});
            % to properly see color differences we need to get rid of extreme values
            areas(areas>mean(areas)+n*std(areas))= mean(areas) + n* std(areas);
            draw_colored_boxes(section_map{zix}, areas, area_bounds, ['Deformation : ' num2str(unique_z(zix))]); % generate figure for y residuals
        end
        
        if options.show_residuals
            % to properly see color differences we need to get rid of extreme values
            tres = tile_residuals_vector{zix};  % all tile residuals for section zu1(zix)
            c = zeros(numel(tres),1);
            for tix = 1:numel(tres)
                c(tix) = sum(tres{tix}(:));
            end
            c(c>mean(c)+n*std(c))= mean(c) + n* std(c);
            draw_colored_boxes(section_map{zix}, c, [0 max(c)], ['Residuals for: ' num2str(unique_z(zix))]); % generate figure for y residuals
            
            %              resx = Resx_store{zix};
            %             resx(resx>mean(resx)+n*std(resx))= mean(resx) + n* std(resx);
            %              resy = Resy_store{zix};
            %             resy(resy>mean(resy)+n*std(resy))= mean(resy) + n* std(resy);
            %
            %resx(resx>maxres) = maxres;
            % resy(resy>maxres) = maxres;
            %resx_bounds(2) = maxres;
            %resy_bounds(2) = maxres;
            %           draw_colored_boxes(sctn_map{zix}, resx, resx_bounds, ['Residuals x: ' num2str(zu1(zix))]); % generate figure for y residuals
            %draw_colored_boxes(sctn_map{zix}, resy, resy_bounds, ['Residuals y: ' num2str(zu1(zix))]); % generate figure for y residuals
        end
    end
end

%% list section outliers and statistics per section
disp('Area should be close to 1.0');
MeanA = zeros(numel(unique_z),1);
MedianA = zeros(numel(unique_z),1);
MedianResidual = zeros(numel(unique_z),1);
MedianPerimeter = zeros(numel(unique_z),1);
ResidualOutliers  = zeros(numel(unique_z),1);
Outliers          = zeros(numel(unique_z), 1);
AreaOutliers  = zeros(numel(unique_z),1);

for zix = 1:numel(unique_z)  % loop over sections
    areas = abs(tile_areas{zix}./(height{zix}.*width{zix}));
    perim = abs(tile_perimeters{zix});
    
    MeanA(zix) = mean(areas);
    MedianA(zix) = median(areas);
    MedianResidual(zix) = confidence(zix);
    MedianPerimeter(zix) = median(perim);
    ResidualOutliers(zix) = numel(residual_outliers_tile_ids{zix});
    AreaOutliers(zix) = numel(area_ratio_outliers_tile_ids{zix});
    Outliers(zix) = numel(outliers_tile_ids{zix});
    SectionName{zix} = num2str(unique_z(zix));
%     str1 = [num2str(zu1(zix)) ' Mean A: ' num2str(mean(areas)) ...
%         ' ----- Median A: ' num2str(median(areas)) ' ---- Median Residual (px): '...
%         num2str(confidence(zix)) ' Median Perimeter: ' num2str(median(perim))];
%     disp([str1]);
    
end
Table = table(MeanA,MedianA,MedianResidual,MedianPerimeter, ResidualOutliers, AreaOutliers, Outliers,...
    'RowNames',SectionName');
disp(Table);
%% generate tile-based residual measure and potential outliers
% section_conf = {};
% for zix = 1:numel(zu1)  % loop over sections
%     
%     tres = res_tiles_vec{zix};  % all tile residuals for section zu1(zix)
%     c = zeros(numel(tres),1);
%     % for each tile, calculate the mean of means of point-match residuals
%     for tix = 1:numel(tres)
%         c(tix) = sum(tres{tix}(:))/numel(tres{tix}(:));
%     end
%     
%     section_conf{zix} = c;
%     
%     % filter by deviation from mean
%     %c(c>mean(c)+n*std(c))= mean(c) + n* std(c);
%     %% determine outliers
%     meann = mean(section_conf{zix});
%     stdd = std(section_conf{zix});
%     I = bsxfun(@gt, abs(bsxfun(@minus, section_conf{zix}, meann)), 2*stdd);
%     outliers_tid{zix} = tidsvec{zix}(find(I));
%     
% end
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
% for zix = 1:1%numel(zu1)
%     c = confidence{zix};
%     t = tidsvec{zix};
%     for tix = 1:size(c,1)
%         if isempty(c{tix})
%             zconf(tix,:) = [nan nan];
%         else
%             zconf(tix,:) = [sum(c{tix},1)/size(c{tix},1)];
%         end
%         disp([num2str(zix) ' '  num2str(tix) ' ' t{tix} '  ' num2str(zconf(tix,:))]);
%     end
% end

%[obj, h, rh, A, minconf, maxconf] = show_map_color(obj, parm, conf,  minconf, maxconf);

%% look at section map
% c = [];
% Ao = 2560*2160;
% colormap jet;
% for zix = 45:45 %numel(zu1)
%     %figure;
%     sm = sctn_map{zix};
%     areas = tile_areas{zix}/Ao;
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













