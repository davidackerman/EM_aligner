function montage_section_SL_prll(fn)
% Intended for deployment: 
% generate montage based on json input provided by fn


% read json input
sl = loadjson(fileread(fn));

if sl.verbose,
    disp('Section montage process started');
    kk_clock();
    disp(['-------  Using input file: ' fn]);
    disp('-------  Using solver options:');disp(sl.solver_options);
    disp('-------  Using solver options:');disp(sl.SURF_options);
    disp('-------  Using source collection:');disp(sl.source_collection);
    disp('-------  Using target collection:');disp(sl.target_collection);    
    disp(['Montage section:' num2str(sl.section_number)]);
    disp('-------  Using target point-match collection:');disp(sl.target_point_match_collection);
end


%%%%%%% start parallel pool
 delete(gcp('nocreate'));
 if sl.verbose, disp('Starting parallel pool');end
% 
      defaultProfile = parallel.defaultClusterProfile;
      myCluster = parcluster(defaultProfile);
      parpool(myCluster);
%       parpool(12);
       poolobj = gcp('nocreate');
 if sl.verbose, disp(['Parallel pool created. Pool size: ' num2str(poolobj.NumWorkers)]);end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 
 
 
if sl.verbose, disp('----- constructing section object');end
L                = Msection(sl.source_collection, sl.section_number);   % instantiate Msection object using the Renderer service to read tiles
L.dthresh_factor = sl.solver_options.dthresh_factor;%1.2;                             % factor x tile diagonal = search radius. increase this paramter to cover a wider radius of tile-tile comparisons

%%% update tile configuration
for tix = 1:numel(L.tiles)
    L.tiles(tix).dir_temp_render = sl.scratch;
    L.tiles(tix).renderer_client = sl.renderer_client;
    L.tiles(tix).fetch_local = 0;
    L.tiles(tix).SURF_NumOctaves = sl.SURF_options.SURF_NumOctaves;
    L.tiles(tix).SURF_NumScaleLevels = sl.SURF_options.SURF_NumScaleLevels;
    L.tiles(tix).SURF_MetricThreshold = sl.SURF_options.SURF_MetricThreshold;
    L.tiles(tix).SURF_MaxFeatures = sl.SURF_options.SURF_MaxFeatures;
end
if sl.verbose
    disp('section data imported --- registration initiated');
    disp(['Found ' num2str(numel(L.tiles)) ' tiles']);
end
if sl.verbose,
tic;
end
%% register
[mL, js]          = register(L, sl.solver_options);                    % perform the actual registration
if sl.verbose
    disp('Section montage finished');
end
if sl.verbose,
toc
end

%% ingest point-matches into point-match database
pm = sl.target_point_match_collection;
if sl.verbose
    disp('Ingesting point matches');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% ingest js into point matches database
    %%% this needs to be done using webwrite --- sosi ---  until then <sigh> we will use curl
    fn = [sl.scratch '/temp_' num2str(randi(100000)) '_' num2str(L.z) '.json'];
    fid = fopen(fn, 'w');
    fwrite(fid,js);
    fclose(fid);
    urlChar = sprintf('%s/owner/%s/matchCollection/%s/matches/', ...
        pm.server, pm.owner, pm.match_collection);
    cmd = sprintf('curl -X PUT --connect-timeout 30 --header "Content-Type: application/json" --header "Accept: application/json" -d "@%s" "%s"',...
        fn, urlChar);
    [a, resp]= evalc('system(cmd)');

%% ingest into Renderer database (optional);
%delete_renderer_stack(sl.target_collection); 
if ~(stack_exists(sl.target_collection))
    disp('Target stack does not exist: creating and then ingesting');
    create_renderer_stack(sl.target_collection);
end
if sl.verbose
    disp('Ingesting results into collection');
end
ingest_section_into_LOADING_collection(mL, sl.target_collection,...
                                       sl.source_collection, pwd, 1); % ingest
resp = set_renderer_stack_state_complete(sl.target_collection);  % set to state COMPLETE

if sl.verbose
    %disp(resp);
    disp('Finished montage and ingestion');
    kk_clock();
end




%%
if sl.verbose, disp('Deleting parallel pool');end
delete(poolobj);

















