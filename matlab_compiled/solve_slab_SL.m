function solve_slab_SL(fn)
% Intended for deployment: solve matrix system based on json input provided by fn


% read json input
sl = loadjson(fileread(fn));

if sl.verbose,
    kk_clock();
    disp(['Using input file: ' fn]);
    disp(['First section:' num2str(sl.first_section)]);
    disp(['Last section:' num2str(sl.last_section)]);
    disp('Using solver options:');disp(sl.solver_options);
    disp('Using source collection:');disp(sl.source_collection);
    disp('Using target collection:');disp(sl.target_collection);
end


[mL, pm_mx, err, R, ~, ntiles, PM, sectionId_load, z_load] = ...
    solve_slab(sl.source_collection, sl.source_point_match_collection, ...
    sl.first_section, sl.last_section, [], sl.solver_options);

%% ingest into Renderer database (optional);
delete_renderer_stack(sl.target_collection);  % delete existing collection if present
ingest_section_into_LOADING_collection(mL, sl.target_collection,...
                                       sl.source_collection, pwd, 1); % ingest
resp = set_renderer_stack_state_complete(sl.target_collection);  % set to state COMPLETE

if sl.verbose
    disp(resp);
    disp('Finished:');
    kk_clock();
end