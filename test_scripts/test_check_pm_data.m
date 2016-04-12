%%%% check point-match data in a collection

%%% check dmesh point-match collection for FAFBv12
clc;kk_clock;
%% configurations

% nfirst = 11;
% nlast  = 16;
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


% nfirst = 1;
% nlast  = 4;
% % configure rough collection
% rc.stack          = ['EXP_v12_rough_' num2str(nfirst) '_' num2str(nlast)];
% rc.owner          ='flyTEM';
% rc.project        = 'test';
% rc.service_host   = '10.37.5.60:8080';
% rc.baseURL        = ['http://' rc.service_host '/render-ws/v1'];
% rc.verbose        = 1;
% 
% % configure point-match collection
% pm.server           = 'http://10.40.3.162:8080/render-ws/v1';
% pm.owner            = 'flyTEM';
% pm.match_collection = 'FAFBv12Test18';


% nfirst = 1;
% nlast  = 16;
% % configure rough collection
% rc.stack          = ['EXP_v12_alignP1_' num2str(nfirst) '_' num2str(nlast)];
% rc.owner          ='flyTEM';
% rc.project        = 'test';
% rc.service_host   = '10.37.5.60:8080';
% rc.baseURL        = ['http://' rc.service_host '/render-ws/v1'];
% rc.verbose        = 1;
% 
% % configure point-match collection
% pm.server           = 'http://10.40.3.162:8080/render-ws/v1';
% pm.owner            = 'flyTEM';
% pm.match_collection = 'FAFBv12Test18';

%%
nfirst = 1;
nlast  = 100;
% configure collection
rc.stack          = ['EXP_v12_rough_' num2str(nfirst) '_' num2str(nlast)];
rc.owner          ='flyTEM';
rc.project        = 'test';
rc.service_host   = '10.37.5.60:8080';
rc.baseURL        = ['http://' rc.service_host '/render-ws/v1'];
rc.verbose        = 1;

% configure point-match collection
pm.server           = 'http://10.40.3.162:8080/render-ws/v1';
pm.owner            = 'flyTEM';
pm.match_collection = 'FAFBv12_set_01';

% configure point match fetching
opts.min_points = 10;
opts.nbrs = 3;
opts.xs_weight = 1/100;
opts.stvec_flag = 1;   % i.e. do not assume rcsource providing the starting values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[L, tIds, PM, pm_mx, sectionIds, zvals] = check_pm_data(nfirst, nlast, rc, pm, opts);