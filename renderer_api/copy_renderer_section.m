function resp = copy_renderer_section(z, rcfrom, rcto, dir_temp)
% copies all tiles with z-value(s) z from rcfrom to rcto
% z can be more than one section

% check
check_input(rcto);
check_input(rcfrom);
if ~stack_exists(rcto)
    disp('Target collection not found, creating new collection in state: ''Loading''');
    resp = create_renderer_stack(rc_target);
end

% configure
verbose = 1;
translate_to_positive_space = 0;
complete = 0;
disableValidation = 0;

parfor zix = 1:numel(z)
    if verbose, disp(['copying section ' num2str(z(zix)) ' from ' rcfrom.stack ' to ' rcto.stack]);end
    
    L = Msection(rcfrom, z(zix)); % read the section with z-value z
    
    delete_renderer_section(rcto, z(zix), 0);
    
    ingest_section_into_renderer_database(L, rcto, rcfrom, dir_temp, ...
        translate_to_positive_space, complete, disableValidation);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_input(rc)
if ~isfield(rc, 'baseURL'), disp_usage; error('baseURL not provided');end
if ~isfield(rc, 'owner'), disp_usage; error('owner not provided');end
if ~isfield(rc, 'project'), disp_usage; error('project not provided');end
if ~isfield(rc, 'stack'), disp_usage; error('stack not provided');end


%%
function disp_usage()
disp('Usage:');
disp('Provide an input struct with fields: baseURL, owner, project, stack');