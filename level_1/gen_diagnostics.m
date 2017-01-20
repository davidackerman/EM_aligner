function [mA, mS, sctn_map, confidence, tile_areas, tile_perimeters, tidsvec, Resx,Resy] =...
    gen_diagnostics(rc, zstart, zend, pm, opts)
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
    pm.server           = 'http://10.40.3.162:8080/render-ws/v1';
    pm.owner            = 'flyTEM';
    pm.match_collection = 'v12_dmesh';
    
end
if nargin<5
    opts.min_points = 5;
    opts.nbrs = 4;
    opts.show_deformation = 1;
    opts.show_residuals = 1;
end


[zu1, sID1, sectionId1, z1, ns1] = get_section_ids(rc, zstart, zend);

%% find corresponding tile centers
tic
% initialize variables to store deformation and residuals
confidence  = cell(numel(zu1),1);
sctn_map  = cell(numel(zu1),1);
tile_areas = cell(numel(zu1),1);
tile_perimeters = cell(numel(zu1),1);
tidsvec = cell(numel(zu1),1);
Resx = cell(numel(zu1),1);
Resy = cell(numel(zu1),1);
res_tiles_vec = cell(numel(zu1),1);
confidence = zeros(numel(zu1,1));
H = cell(numel(zu1),1);
W = cell(numel(zu1),1);
% to generate histogram counts we need to define bin edges
edges = [0.4:.02:1.7];
counts = zeros(numel(zu1), numel(edges));
webopts = weboptions('Timeout', 60);
for zix = 1:numel(zu1)
    Ar = [];    % surface area of tiles
    Aratio = [];% surface area ratio of tiles
    S  = [];    % perimeter of tiles
    p = [];     % tile patch
    height = [];
    width = [];
    % call the Renderer API to fetch tile information
    urlChar = sprintf('%s/owner/%s/project/%s/stack/%s/z/%.1f/tile-specs', ...
        rc.baseURL, rc.owner, rc.project, rc.stack,zu1(zix) );
    j = webread(urlChar, webopts);
    
    % AREA and PERIMETER process individual tiles: Calculate area and perimeter
    jt1 = tile;
    sectionID = j(1).layout.sectionId;
    count = 1;
    tids = {};
    for jix = 1:numel(j)
        jt1(jix) = tile(j(jix));
        tids(jix) = {jt1(jix).renderer_id};
        % make four corners for this tile
        x = 0;
        y = 0;
        Px = [x; x + jt1(jix).W; x + jt1(jix).W; x];
        Py = [y; y    ; y + jt1(jix).H; jt1(jix).H];
        %%% transform points
        if strcmp(class(jt1(jix).tform), 'affine2d')
            P = [Px(:) Py(:) [1 1 1 1]']*jt1(jix).tform.T;
        else
            P = transformPointsInverse(jt1(jix).tform,[Px Py]);
        end
        p{jix} = {P};
        %cm(jix,:) = sum([P(:,1)/4 P(:,2)/4],1);  % center of mass of each tile
        
        height(count) = jt1(jix).H;
        width(count)  = jt1(jix).W;
        % check polygon area
        Ar(count) = polyarea(P(:,1), P(:,2));
        Aratio(count) = Ar(count)/(jt1(jix).H * jt1(jix).W);
        
        %%% polygonperimeter
        s = 0;
        s = s + sqrt((P(1,1)-P(2,1)).^2 + (P(1,2)-P(2,2)).^2);
        s = s + sqrt((P(2,1)-P(3,1)).^2 + (P(2,2)-P(3,2)).^2);
        s = s + sqrt((P(3,1)-P(4,1)).^2 + (P(3,2)-P(4,2)).^2);
        s = s + sqrt((P(1,1)-P(4,1)).^2 + (P(1,2)-P(4,2)).^2);
        S(count) = s;
        
        
        count = count + 1;
    end
    H{zix} = height;
    W{zix} = width;
    tidsvec{zix} = tids;
    counts(zix,:) = histc(Aratio, edges, 2);
    tile_areas{zix} = Ar;
    tile_perimeters{zix} = S;
    mA(zix) = median(Aratio);
    mS(zix) = median(S);
    sctn_map{zix} = p;   % needed to plot tile boxes
    
    %% %%% determine point-matches, solution and residuals for this section
    
    % First: load point-matches and section into "L" (point-matches are in L's pm struct field)
    if (zix + opts.nbrs+1)>numel(zu1)
        pmz2 = zu1(end);
    else
        pmz2 = zu1(zix + opts.nbrs);
    end
    [L]  = ...
        load_point_matches(zu1(zix), pmz2, rc, pm, opts.nbrs, ...
        opts.min_points, 0);
    
    % initialize variable tpr to hold point-match residuals
    tpr = zeros(numel(L.tiles),1);   % tpr === transformed point-match residuals
    res_tiles = cell(numel(L.tiles),1);
    % Second: generate point-match residuals from L.pm by transforming them and taking the sum of squared
    % residuals
    res_vecx = [];
    res_vecy = [];
    for pmix = 1:size(L.pm.M,1)
        a1 = L.pm.adj(pmix,1);
        a2 = L.pm.adj(pmix,2);
        m1 = L.pm.M{pmix,1};
        m2 = L.pm.M{pmix,2};
        m1 = [m1 ones(size(m1,1),1)]*L.tiles(a1).tform.T;  % apply transformation
        m2 = [m2 ones(size(m2,1),1)]*L.tiles(a2).tform.T;  % apply transformation
        res = real(sqrt(sum((m1(1)-m2(1))*(m1(1)-m2(1))  + (m1(2)-m2(2)) * (m1(2)-m2(2)))));    %%%% sum of squared residuals
        %disp(res);
        %res_vec(pmix,:) = sqrt((m1(1)-m2(1))*(m1(1)-m2(1)) + (m1(1)-m2(1))*(m1(2)-m2(2)));
        
        %%%% sosi
        %         res_vecx(pmix,:) = abs(m1(1)-m2(1));
        %         res_vecy(pmix,:) = abs(m1(2)-m2(2));
        %
        tpr(a1) = [tpr(a1) + res];
        tpr(a2) = [tpr(a2) + res];
        res_tiles{a1} = [res_tiles{a1} res];
        res_tiles{a2} = [res_tiles{a2} res];
    end
    res_tiles_vec{zix} = res_tiles;  % store tile residuals for this section
    Resx{zix} = res_vecx;
    Resy{zix} = res_vecy;
    % store errors/confidence in cell array "confidence"
    confidence(zix) = 0;
    for tix = 1:numel(L.tiles)
        confidence(zix) = confidence(zix) + tpr(tix);%{sum(tpr{tix},1)/size(tpr{tix},1)};
    end
    
end
Resx_store = Resx;
Resy_store = Resy;
Resx = cell2mat(Resx);
Resy = cell2mat(Resy);
%%


%%
if opts.show_deformation || opts.show_residuals
    figure;plot(zu1, confidence, '-o', 'LineWidth',2);title(['Point-match residuals: ' num2str(zstart) ' to ' num2str(zend)]);
    xlabel('Section number');
    ylabel('Point-match residual (sum of sum of squared distance)');
    figure;plot(zu1, mA, '-o', 'LineWidth',2);title(['Area ratio median per layer: ' num2str(zstart) ' to ' num2str(zend)]);
    xlabel('Section number');
    ylabel('Area ratio median (ideally one)');
    %figure;plot(zu1, mS, 'LineWidth',2);title('Perimeter median perlayer');
    n = 0.1; % number of std for cutoff to determine outliers
    maxres = 2.0;  % cap res at this value so that color differences can be discerned
    
    Resx(Resx>mean(Resx)+n*std(Resx))= mean(Resx) + n* std(Resx);
    Resy(Resy>mean(Resy)+n*std(Resy))= mean(Resy) + n* std(Resy);
    
    
    % figure; hist(real(Resx), 100); title([char(rc.stack) ' -- histogram for all tiles Resx']);axis tight
    % figure; hist(real(Resy), 100); title([char(rc.stack) ' -- histogram for all Resy']);axis tight
    
    %% display section tile box figures
    % at this point Resx and Resy contain residual information
    % tile_areas contains areas for individual tiles
    areas = [cell2mat(tile_areas')]'./[cell2mat(H').*cell2mat(W')]';
    areas = abs(1-areas);  % measure area deviation from 1
    %%% display histogram of tile areas for full slab
    figure;hist(areas,100);title('Deformation histogram: (ideally zero) for whole slab');axis tight
    xlim([0 1]);
    c_bounds = [min(0) max(20)];
    area_bounds = [0 2];
    colormap jet;
    
    for zix = 1:numel(zu1)  % loop over sections
        if opts.show_deformation
            areas = tile_areas{zix}./(H{zix}.*W{zix});
            % to properly see color differences we need to get rid of extreme values
            areas(areas>mean(areas)+n*std(areas))= mean(areas) + n* std(areas);
            draw_colored_boxes(sctn_map{zix}, areas, area_bounds, ['Deformation : ' num2str(zu1(zix))]); % generate figure for y residuals
        end
        
        if opts.show_residuals
            % to properly see color differences we need to get rid of extreme values
            tres = res_tiles_vec{zix};  % all tile residuals for section zu1(zix)
            c = zeros(numel(tres),1);
            for tix = 1:numel(tres)
                c(tix) = sum(tres{tix}(:));
            end
            c(c>mean(c)+n*std(c))= mean(c) + n* std(c);
            draw_colored_boxes(sctn_map{zix}, c, [0 max(c)], ['Residuals for: ' num2str(zu1(zix))]); % generate figure for y residuals
            
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
for zix = 1:numel(zu1)  % loop over sections
    areas = abs(tile_areas{zix}./(H{zix}.*W{zix}));
    disp('Area should be close to 1.0');
    str1 = [num2str(zu1(zix)) ' dA -Mean: ' num2str(mean(areas)) ...
        ' ----- Median: ' num2str(median(areas)) ' ---- Residual: ' num2str(confidence(zix))];
    %     resx = Resx_store{zix};
    %     resy = Resy_store{zix};
    %
    %     str2 = [ ' resx: Median: ' num2str(median(resx)) ' -outliers' num2str(numel(find(resx>maxres)))];
    %     str3 = [ ' resy: Median: ' num2str(median(resy)) ' -outliers' num2str(numel(find(resy>maxres)))];
    disp([str1]);
    
    %     resx(resx>maxres) = maxres;
    %     resy(resy>maxres) = maxres;
    
    %     figure;
    %     subplot(3,1,1); hist(areas, 100); title([num2str(zu1(zix)) '-- Area deviation']); xlim([0 0.1]);
    %     subplot(3,1,2);  hist(resx, 100); title([num2str(zu1(zix)) '-- resx']); xlim([0 maxres]);
    %     subplot(3,1,3); hist(resy, 100); title([num2str(zu1(zix)) '-- resy']); xlim([0 maxres]);
    %     drawnow;
end

%% summarize deformation for whole stack
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













