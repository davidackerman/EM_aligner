%% This is a test script for generation of the montage of one section
% The steps are:
% [0] configuration of renderer database collection (source and output montage)
% [1] generate the Msection object and register (generate montage)
% [2] ingest the montaged section into the Renderer database as a new collection
% [3] quick inspection of results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% [0] configure collections

section_z = 4.0; % this is the section z-coordinate value to be montaged


% configure source collection
rcsource.stack          = 'v12_acquire_merged';
rcsource.owner          = 'flyTEM';
rcsource.project        = 'FAFB00';
rcsource.server         = 'http://10.37.5.60:8080/render-ws/v1'; % use of ip adress is preferred (no DNS lookup)-- Note: 10.37.5.60 is a VM, 10.40.3.162 is tem-services
rcsource.service_host   = '10.37.5.60:8080';
rcsource.verbose        = 1;
rcsource.baseURL        = rcsource.server;
rcsource.source_stack   = rcsource.stack;
rcsource.verbose        = 1;

% configure output montage target collection
rctarget_montage = rcsource;            % initialize to the same as above
rctarget_montage.project = 'test';      % we will put the montage under "test" projects
rctarget_montage.collection = ['EXP_' rcsource.stack '_' num2str(section_z)];


%% [1] generate montage registration
L                = Msection(rcsource, section_z);   % instatiate an Msection object using the Renderer service to read tiles
L.dthresh_factor = 0.9;                             % factor x tile diagonal = search radius. increase this paramter to cover a wider radius of tile-tile comparisons
[L, ~]           = register(L);                     % perform the actual registration

%% [2] ingest montaged section into the database
ingest_section_into_renderer_database(mL, rctarget_montage, rcsource, collection, pwd);
mL = update_tile_sources(mL, rctarget_montage.owner, rctarget_montage.project, collection, rcsource.server); % point to the new collection

%% [3] look at result in Matlab
figure; show_map(L);title('Tile outlines before registration');
figure; show_map(mL);title('Tile outlines after registration');


