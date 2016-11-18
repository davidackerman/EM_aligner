%% test jdr-based montage solve
fn = '/groups/flyTEM/home/khairyk/EM_aligner/test_scripts/test_jdr_solve_montage_input.json';
% read json input
sl = loadjson(fileread(fn));

if sl.verbose
    kk_clock();
    disp(sl);
    disp(['Using input file: ' fn]);
    disp(['Section(s) with z value:' num2str(sl.z_value)]);
    disp('Using solver options:');disp(sl.solver_options);
    disp('Using jdr-specific solver options:');disp(sl.solver_options.jdr_options);
    disp('Using source collection:');disp(sl.source_collection);
    disp('Using target collection:');disp(sl.target_collection);
end
%%
tic;if sl.verbose, disp('-- Loading point matches');end
[L, tIds, PM, pm_mx, sectionId_load, z_load]  = ...
    load_point_matches(sl.z_value,sl.z_value, sl.source_collection, ...
    sl.source_point_match_collection, 0, sl.solver_options.min_points, 0, sl.solver_options.max_points); %
toc
if sl.filter_point_matches
    tic;if sl.verbose, disp('-- Filtering point matches');end
    pmfopts.NumRandomSamplingsMethod = 'Desired confidence';
    pmfopts.MaximumRandomSamples = 1000;
    pmfopts.DesiredConfidence = 99.9;
    pmfopts.PixelDistanceThreshold = 0.1;
    if sl.verbose,
        disp('using point-match filter:');
        disp(pmfopts);
    end
    L.pm = filter_pm(L.pm, pmfopts);
    toc
end
%% add translation peggs for matrix stability
if sl.solver_options.use_peg
    L = add_translation_peggs(L, sl.solver_options.peg_npoints, sl.solver_options.peg_weight);
end
[L, ntiles] = reduce_to_connected_components(L, sl.solver_options.min_tiles);
L = L(1);
%% solve rigid approximation
disp('Generating json rotation approximation tile-spec file ...');
sl.solver_options.distributed = 0;
[Lr, errR, mL, is, it, Res, ...
    A, b, B, d, W, K, Lm, xout, L2, U2, tB, td,...
    At, bt, Bt, dt, Wt, Kt, Lmt, xoutt, L2t, U2t, tBt, tdt]  =...
    get_rigid_approximation(L, 'backslash', sl.solver_options);

if sl.solver_options.use_peg
    %%% remove peggs and last tile
    last_tile = numel(Lr.tiles);
    del_ix = find(Lr.pm.adj(:,2)==last_tile);
    Lr.pm.M(del_ix,:)  = [];
    Lr.pm.adj(del_ix,:) = [];
    Lr.pm.W(del_ix) = [];
    Lr.pm.np(del_ix) = [];
    Lr.tiles(end) = [];
    Lr = update_adjacency(Lr);
end
%% sosi ---- solve affine and export expected transofmormations to compare

[mL, err1, Res1, A, b, B, d, W, K, Lm, xout, LL2, U2, tB, td,...
    invalid, time_Axb, time_gen_A] = ...
   solve_affine_explicit_region(Lr, sl.solver_options);

dir_scratch = '/groups/flyTEM/home/khairyk/solver_paper_work/scratch';
fn_canvas_affine_backslash = [dir_scratch '/tmp_affine_backslash_canvases.json'];
jstr_aff_backslash = export_json(mL, fn_canvas_affine_backslash);
disp(jstr_aff_backslash);
[mLc, tpr, resout_backslash] = tile_based_point_pair_errors(mL, A, xout);
mean_error_affine_backslash = mean(resout_backslash);
if sl.solver_options.use_peg
    %%% remove peggs and last tile
    last_tile = numel(mL.tiles);
    del_ix = find(mL.pm.adj(:,2)==last_tile);
    mL.pm.M(del_ix,:)  = [];
    mL.pm.adj(del_ix,:) = [];
    mL.pm.W(del_ix) = [];
    mL.pm.np(del_ix) = [];
    mL.tiles(end) = [];
    mL = update_adjacency(mL);
end
%%
dir_scratch = '/groups/flyTEM/home/khairyk/solver_paper_work/scratch';
delete([dir_scratch '/*.json']);
disp('calling jdr solver');
debugjdr = 1;
lambda = sl.solver_options.lambda;
degree = sl.solver_options.degree;

fnpmjson = ...
    [dir_scratch '/tmp_pm.json'];
disp('Generating json point-match file ...');
tic;jstr = PM_json(L, fnpmjson);toc;
fn_canvas_json_output = ...
    [dir_scratch '/canvases_out.json'];
solv_cmd = '/groups/flyTEM/home/khairyk/downloads/jdr/joint-image-registration-solver/jdr_solver';

if sl.solver_options.jdr_options.stvec==0  % calculate rigid_approximation followed by affine
    % works!
    disp('JDR: calculating both rigid and affine');
    L_jdr = L;
    stvec = 0;
    fn_canvas_input = '';
else                                      % calcualte affine based on pre-calculated rigid
    % DOES NOT WORK YET --- DEBUG PLEASE
    disp('Calculating only affine based on provided rigid approximation');
    L_jdr = Lr;
    stvec = 1;   % use fn_canvas_input as input for regularization
    fn_canvas_input = [dir_scratch '/tmp_rap_canvases.json'];
    jstr_rap = export_json(L_jdr, fn_canvas_input);
    if debugjdr
        type(fn_canvas_input);
    end
end

time_Axb = -999;
time_gen_A = -999;
%%
disp('-------------------- Invoking cpp solver --------------');
cmd = [solv_cmd ' ' fnpmjson ' ' fn_canvas_json_output ' ' num2str(degree) ...
    ' ' num2str(lambda) ' ' num2str(stvec) ' ' fn_canvas_input];
tic;[a,resp_str] = system(cmd);toc

%%% diagnostic and profiling information

if debugjdr
    disp(resp_str);
    disp(fn_canvas_json_output);
    type(fn_canvas_json_output);
    disp('expected:');
    disp(jstr_aff_backslash);
end

if degree>0
C = strsplit(resp_str, 'TIME_');
c = strsplit(C{2}, ' ');
time_gen_A = str2double(c{2});
c = strsplit(C{3}, ' ');
time_Axb = str2double(c{2});
end
disp('finished cpp solution');
kk_clock;
%disp(resp_str);
if sl.solver_options.jdr_options.stvec==0
    mLjdr = update_transformation_from_json(L,fn_canvas_json_output);
    
    %%% remove peggs and last tile
    last_tile = numel(mLjdr.tiles);
    del_ix = find(mLjdr.pm.adj(:,2)==last_tile);
    mLjdr.pm.M(del_ix,:)  = [];
    mLjdr.pm.adj(del_ix,:) = [];
    mLjdr.pm.W(del_ix) = [];
    mLjdr.pm.np(del_ix) = [];
    mLjdr.tiles(end) = [];
    mLjdr = update_adjacency(mLjdr);
else
    
    mLjdr = update_transformation_from_json(Lr,fn_canvas_json_output);
end
%% % to continue diagnostics we need to set to matrix_only
sld = sl;
sld.solver_options.matrix_only = 1;
[mL2, err2, Res1, Ajdr, b, B, d, W, K, Lm, xoutjdr, LL2, U2, tB, td,...
    invalid] = solve_affine_explicit_region(mLjdr,...
    sld.solver_options);
[mLc, tpr, resout_jdr] = tile_based_point_pair_errors(mLjdr, Ajdr, xoutjdr);
mean_error_affine_jdr = mean(resout_jdr);

%%
%% ingest into Renderer database (optional);
if sl.target_collection.initialize,
    if sl.verbose, disp('Initializing collection / Deleting existing');end
    delete_renderer_stack(sl.target_collection);  % delete existing collection if present
end
tic;if sl.verbose, disp('-- Ingesting section into collection');end
resp = ingest_section_into_LOADING_collection(mLjdr, sl.target_collection,...
    sl.source_collection, sl.temp_dir, 1); % ingest
if sl.target_collection.complete
    if sl.verbose, disp('Completing collection');end
    resp = set_renderer_stack_state_complete(sl.target_collection);  % set to state COMPLETE
end
if sl.verbose
    disp(resp);
    disp('Finished:');
    kk_clock();
end

%% ingest reference backslash solution
rc = sl.target_collection;
rc.stack = 'EXP_test_montage_solver_backslash';
if sl.target_collection.initialize
    if sl.verbose, disp('Initializing collection / Deleting existing');end
    delete_renderer_stack(rc);  % delete existing collection if present
end
tic;if sl.verbose, disp('-- Ingesting section into collection');end
resp = ingest_section_into_LOADING_collection(mL, rc,...
    sl.source_collection, sl.temp_dir, 1); % ingest
if sl.target_collection.complete
    if sl.verbose, disp('Completing collection');end
    resp = set_renderer_stack_state_complete(rc);  % set to state COMPLETE
end
if sl.verbose
    disp(resp);
    disp('Finished:');
    kk_clock();
end
%% optional
str = view_collection_dashboard(sl.target_collection); disp(str);

%%
clc
if debugjdr
    
    
    fn_ = [dir_scratch '/tmp_rap_canvases.json'];
    jstr_rap = export_json(Lr, fn_);
    disp('Rigid approximation');
    type(fn_canvas_input);
    disp([ ' JDR: output: ' fn_canvas_json_output]);
    type(fn_canvas_json_output);
    disp('expected:');
    disp(jstr_aff_backslash);
    disp([err1 err2]);
    disp([mean_error_affine_backslash mean_error_affine_jdr]);
end
gen_section_based_tile_deformation_statistics(sl.target_collection,...
    sl.z_value, sl.z_value, sl.source_point_match_collection);

gen_section_based_tile_deformation_statistics(rc,...
    sl.z_value, sl.z_value, sl.source_point_match_collection);
