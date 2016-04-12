function [L, tIds, PM, pm_mx, sectionIds, zvals]  = check_pm_data(nfirst, nlast, rc, pm, opts)
%%%% check point-match data in a collection visually
% example usage:
%%% check dmesh point-match collection for FAFBv12
% clc;kk_clock;
% nfirst = 100;
% nlast  = 105;
%
%
% % configure rough collection
% rc.stack          = 'v12_align';
% rc.owner          ='flyTEM';
% rc.project        = 'FAFB00';
% rc.service_host   = '10.37.5.60:8080';
% rc.baseURL        = ['http://' rc.service_host '/render-ws/v1'];
% rc.verbose        = 1;
%
% % configure point-match collection
% pm.server           = 'http://10.40.3.162:8080/render-ws/v1';
% pm.owner            = 'flyTEM';
% pm.match_collection = 'v12_dmesh';
%
% % % configure rough collection
% % rc.stack          = ['EXP_v12_rough_' num2str(nfirst) '_' num2str(nlast)];
% % rc.owner          ='flyTEM';
% % rc.project        = 'test';
% % rc.service_host   = '10.37.5.60:8080';
% % rc.baseURL        = ['http://' rc.service_host '/render-ws/v1'];
% % rc.verbose        = 1;
% %
% % % configure point-match collection
% % pm.server           = 'http://10.40.3.162:8080/render-ws/v1';
% % pm.owner            = 'flyTEM';
% % pm.match_collection = 'FAFBv12Test17';
%
%
% opts.min_points = 10;
% opts.nbrs = 3;
% opts.xs_weight = 1/100;
% opts.stvec_flag = 1;   % i.e. do not assume rcsource providing the starting values.
%
%
% Author: Khaled Khairy. Janelia Research Campus
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[L, tIds, PM, pm_mx, sectionIds, zvals] = load_point_matches(nfirst, nlast, rc, pm, opts.nbrs, opts.min_points, 0.1);
L = update_tile_sources(L, rc);
T = array2table(pm_mx);disp(T);

%%% plot cartoons based on rough collection

% associate each tile with a coordinate
for tix = 1:numel(L.tiles)
    z(tix) = L.tiles(tix).z;
    j = get_tile_spec_renderer(L.tiles(tix));
    str = j.tileSpecs.transforms.specList(end).dataString;
    c = strsplit(str);
    x(tix) = str2double(c{5});
    y(tix) = str2double(c{6});
end

% loop over edges of L.G to plot edges
figure(1);
ax = gca;
for eix = 1:size(L.G.Edges,1)
    edg = table2cell(L.G.Edges(eix,1));
    rid1 = edg{1}{1};
    rid2 = edg{1}{2};
    indx1 = L.map_renderer_id(rid1);
    indx2 = L.map_renderer_id(rid2);
    
    X = [x(indx1) x(indx2)];
    Y = [y(indx1) y(indx2)];
    Z = [z(indx1) z(indx2)];
    if Z(1)==Z(2)
        line(X,Y,Z,'LineWidth',1, 'Color',[.5 .5 .5],'Parent',ax);hold on;
    else
        line(X,Y,Z,'LineWidth',1, 'Color',[.2 .2 .2],'Parent',ax);hold on;
    end
    plot3(X,Y,Z, '*b');
end

% loop over sections and generate rectangles
for six = 1:numel(sectionIds)
    [Wbox, bbox, url] = get_section_bounds_renderer(rc, zvals(six));
    X = [bbox(1) bbox(3) bbox(3) bbox(1)];
    Y = [bbox(2) bbox(2) bbox(4) bbox(4)];
    Z = [zvals(six) zvals(six) zvals(six) zvals(six)];
    figure(1);hold on; patch(X, Y, Z, 'red', 'FaceAlpha', 0.2);
end

daspect('auto');
view(3);cameramenu;



























