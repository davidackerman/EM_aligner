function system_solve_translation_SL(fn)
% Intended for deployment: solve matrix system using translation based on json input provided by fn

% read json input
sl = loadjson(fileread(fn));
if ~isfield(sl.solver_options, 'diagnostics_level'), sl.solver_options.diagnostics_level = 0; end
if sl.verbose
    kk_clock();
    disp(['Using input file: ' fn]);
    disp(['First section:' num2str(sl.first_section)]);
    disp(['Last section:' num2str(sl.last_section)]);
    disp('Using solver options:');disp(sl.solver_options);
    disp('Using source collection:');disp(sl.source_collection);
    disp('Using target collection:');disp(sl.target_collection);
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

%%% deprecated?
[err,R, Tout, PM, Diagnostics] = system_solve_translation(sl.first_section, sl.last_section, sl.source_collection, sl.source_point_match_collection, sl.solver_options, sl.target_collection);
if sl.solver_options.diagnostics_level>=0
    fprintf('root_mean_square_residual_value_mean: %f\n',mean(Diagnostics.rms));
    if sl.solver_options.diagnostics_level == 1
        fprintf('root_mean_square_residual_values: ');
        fprintf('%f ',Diagnostics.rms);
        fprintf('\n');
    end
end