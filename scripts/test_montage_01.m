%% The full stitching pipeline working in conjunction with the Renderer database and the point-match database
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [0] configure collections

% configure source collection
rcsource.stack = 'v12_acquire_merged';
rcsource.owner='flyTEM';
rcsource.project= 'FAFB00';
rcsource.server = 'http://10.37.5.60:8080/render-ws/v1';
rcsource.service_host = '10.37.5.60:8080';
rcsource.verbose = 1;
rcsource.baseURL        = rcsource.server;
rcsource.source_stack   = rcsource.stack;
rcsource.verbose        = 1;

% configure montage target collection
rctarget_montage = rcsource;  % initialize to the same as above
rctarget_montage.project = 'test'; % we will put the montage under "test" projects

section_z 4;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% get the list of zvalues and section ids within the z range between nfirst and nlast (inclusive)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
urlChar = sprintf('%s/owner/%s/project/%s/stack/%s/sectionData', ...
    rcsource.server, rcsource.owner, rcsource.project, rcsource.stack);
j = webread(urlChar);
sectionId = {j(:).sectionId};
z         = [j(:).z];
indx = find(z = section_z);

sectionId = sectionId(indx);% determine the sectionId list we will work with
z         = z(indx);        % determine the zvalues (this is also the spatial order)

%% generate montage registration
L = Msection(rcsource, z);
L.dthresh_factor = 1.0;
L = update_XY(L);
L  = update_adjacency(L);
[L, js] = alignTEM_inlayer(L);

%% ingest
opts.outlier_lambda = 1e3;  % large numbers result in fewer tiles excluded
mL = concatenate_tiles(L, opts.outlier_lambda);  % usually only relevant if we have more than one section
% i.e. if L were an Msection array
collection = ['EXP_v12_SURF_montage_' num2str(section_z)];
ingest_section_into_renderer_database(mL, rctarget_montage, rcsource, collection, pwd);
mL = update_tile_sources(mL, rctarget_montage.owner, rctarget_montage.project, collection, rcsource.server); % point to the new collection

