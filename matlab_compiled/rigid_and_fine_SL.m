function rigid_and_fine_SL(fn)
% read json input


sl = loadjson(fileread(fn));
if ~isfield(sl.rigid_solver_options, 'diagnostics_level'), sl.rigid_solver_options.diagnostics_level = 0; end
if ~isfield(sl.fine_solver_options, 'diagnostics_level'), sl.fine_solver_options.diagnostics_level = 0; end

if sl.verbose
    kk_clock();
    disp(['Using input file: ' fn]);
    disp(['First section:' num2str(sl.first_section)]);
    disp(['Last section:' num2str(sl.last_section)]);
    disp('Using rigid solver options:');disp(sl.rigid_solver_options);
    disp('Using fine solver options:');disp(sl.fine_solver_options);
    disp('Using source collection:');disp(sl.source_collection);
    disp('Using rigid collection:');disp(sl.rigid_collection);
    disp('Using fine collection:');disp(sl.fine_collection);
    for i = 1:numel(sl.source_point_match_collection)
        disp('Using point-match collection:');disp(sl.source_point_match_collection(i));
    end
end

if isfield(sl, 'num_cores')
    if sl.num_cores>1
        parpool(sl.num_cores)
    end
else
   parpool 
end
disp('Solve for rigid');
%sl.rigid_collection.versionNotes = gen_versionNotes(sl.rigid_solver_options);
[err,R, Tout, A, b, map_id, tIds, z_val, Diagnostics] = system_solve_rigid_approximation(sl.first_section, sl.last_section, sl.source_collection, sl.source_point_match_collection, sl.rigid_solver_options, sl.rigid_collection);

%%% deprecated?
if sl.rigid_solver_options.diagnostics_level>=0
    fprintf('root_mean_square_residual_value_mean: %f\n',mean(Diagnostics.rms));
    if sl.rigid_solver_options.diagnostics_level == 1
        fprintf('root_mean_square_residual_values: ');
        fprintf('%f ',Diagnostics.rms);
        fprintf('\n');
    end
end

disp('Solve for affine');
sl.fine_collection.versionNotes = gen_versionNotes(sl.fine_solver_options);
[err,R, Tout_affine, D_affine] = system_solve_affine_with_constraint(sl.first_section, sl.last_section, sl.rigid_collection,sl.source_point_match_collection, sl.fine_solver_options, sl.fine_collection);
disp(err);
%%% deprecated?
if sl.fine_solver_options.diagnostics_level>=0
    fprintf('root_mean_square_residual_value_mean: %f\n',mean(D_affine.rms));
    if sl.fine_solver_options.diagnostics_level == 1
        fprintf('root_mean_square_residual_values: ');
        fprintf('%f ',D_affine.rms);
        fprintf('\n');
    end
end
end

