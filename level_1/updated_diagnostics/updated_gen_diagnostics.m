function [output_struct] =...
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
    options.number_of_cross_sections=2;
end

% %%% defaults and overrides
if ~isfield(options, 'show_residual_histogram'), options.show_residual_histogram = 0;end
if ~isfield(options, 'nstd'), options.nstd = 2;end
if ~isfield(options, 'residual_info'), options.residual_info = 0;end

tic;
[ output_struct, all_height, all_width, counts, all_section_map, section_ids_grouped_by_z, unique_z ] = calculate_area_and_perimeter_ratios(rcsource, rc, zstart, zend, options);
toc;
%%
%% Montage and Cross-section residuals
if options.residual_info
    [ output_struct2 ] = calculate_montage_point_match_residuals(rc, point_matches, options, unique_z);
    output_struct.MontageResiduals.residuals = output_struct2.MontageResiduals.residuals;
    output_struct.MontageResiduals.median_of_means = output_struct2.MontageResiduals.median_of_means;
    output_struct.MontageResiduals.number_of_outliers = output_struct2.MontageResiduals.number_of_outliers;
    output_struct.MontageResiduals.outlier_tile_ids = output_struct2.MontageResiduals.outlier_tile_ids;
    output_struct.MontageResiduals.number_of_unconnected_tiles = output_struct2.MontageResiduals.number_of_unconnected_tiles;
    output_struct.MontageResiduals.unconnected_tile_ids = output_struct2.MontageResiduals.unconnected_tile_ids;
    residuals_matrix = calculate_cross_section_point_match_residuals(rc, zstart, zend, point_matches, options,unique_z, section_ids_grouped_by_z);
    residuals_matrix = residuals_matrix + diag(output_struct.MontageResiduals.median_of_means);
    output_struct.CrossSectionAndMontageResidualsMatrix = residuals_matrix;
end
%%
n = 0.1; % number of std for cutoff to determine outliers

%%
if options.show_deformation || options.show_residuals
    figure;plot(unique_z, all_montage_residuals_vector, '-o', 'LineWidth',2);title(['Point-match residuals: ' num2str(zstart) ' to ' num2str(zend)]);
    xlabel('Section number');
    ylabel('Point-match residual (sum of sum of squared distance)');
    figure;plot(unique_z, all_area_ratio_median, '-o', 'LineWidth',2);title(['Area ratio median per layer: ' num2str(zstart) ' to ' num2str(zend)]);
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
    areas = [cell2mat(output_struct.Area.areas')]'./[cell2mat(all_height').*cell2mat(all_width')]';
    areas = abs(1-areas);  % measure area devisation from 1
    %%% display histogram of tile areas for full slab
    figure;hist(areas,100);title('Deformation histogram: (ideally zero) for whole slab');axis tight
    xlim([0 1]);
    c_bounds = [min(0) max(20)];
    area_bounds = [0 2];
    colormap jet;
    
    for z_index = 1:numel(unique_z)  % loop over sections
        if options.show_deformation
            areas = output_struct.Area.areas{z_index}./(all_height{z_index}.*all_width{z_index});
            % to properly see color differences we need to get rid of extreme values
            areas(areas>mean(areas)+n*std(areas))= mean(areas) + n* std(areas);
            draw_colored_boxes(all_section_map{z_index}, areas, area_bounds, ['Deformation : ' num2str(unique_z(z_index))]); % generate figure for y residuals
        end
        
        if options.show_residuals
            % to properly see color differences we need to get rid of extreme values
            tres = output_struct.MontageResiduals.residuals{z_index};  % all tile residuals for section zu1(z_index)
            c = zeros(numel(tres),1);
            for tix = 1:numel(tres)
                c(tix) = sum(tres{tix}(:));
            end
            c(c>mean(c)+n*std(c))= mean(c) + n* std(c);
            draw_colored_boxes(all_section_map{z_index}, c, [0 max(c)], ['Residuals for: ' num2str(unique_z(z_index))]); % generate figure for y residuals
            
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
    MedianAreaRatio = output_struct.Area.median_of_means;
    MedianPerimeterRatio = output_struct.Perimeter.median_of_means;
    MedianMontageResidual = output_struct.MontageResiduals.median_of_means;
    AreaOutliers  = output_struct.Area.number_of_outliers;
    PerimeterOutliers  = output_struct.Perimeter.number_of_outliers;
    MontageResidualOutliers  = output_struct.MontageResiduals.number_of_outliers;
    TotalAreaAndPerimeterOutliers = cellfun(@(area_outlier_tile_ids, perimeter_outlier_tile_ids) numel(unique([area_outlier_tile_ids(:);perimeter_outlier_tile_ids(:)])), output_struct.Area.outlier_tile_ids,  output_struct.Perimeter.outlier_tile_ids);
    SectionName = cellstr(num2str(unique_z'));
    Table = table(MedianAreaRatio,MedianPerimeterRatio,MedianMontageResidual, AreaOutliers, PerimeterOutliers, MontageResidualOutliers, TotalAreaAndPerimeterOutliers,...
        'RowNames',SectionName');
    output_struct.Table = Table;
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












