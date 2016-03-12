function ingest_section_into_renderer_database(mL,rctarget, rcsource, target_stack, dir_work)


%% translate to origin to be Renderer friendly
disp('translating to set in +ve space');
mL = translate_to_origin(mL);
%% export to MET (in preparation to be ingested into the Renderer database
disp('-------------------- Exporting to MET -------------');
zix = mL.z;
fn_MET = [dir_work '/X_A_' num2str(zix) '.txt'];
if strcmp(class(mL.tiles(1).tform), 'images.geotrans.PolynomialTransformation2D')
    export_montage_MET2(mL, fn_MET);
    v = 'v3';
else
    export_MET(mL, fn_MET, 2, 2, 0);
    v = 'v1';
end
%% ingest into the Renderer database
% generate a new collection in the Renderer database


rc = rctarget;

rc.target_project = rctarget.project;
rc.source_project = rcsource.project;
rc.target_stack   = target_stack;
rc.source_stack   = rcsource.stack;

rc.MET_format     = v;
rc.MET_files{1}   = fn_MET;
rc.baseURL = rc.server;
generate_results_collection(rc);