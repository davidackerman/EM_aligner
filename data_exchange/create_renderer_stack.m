function [err, resp] = create_renderer_stack(rc)
% Creates a new collection in the Renderer database
% rc is a struct with fields (baseURL, owner, project, stack)
%
% Author: Khaled Khairy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbose = 0;
err = 0;
check_input(rc);

str1 = sprintf('PROJECT_PARAMS="--baseDataUrl %s --owner %s --project %s";', rc.baseURL, rc.owner, rc.project);
str3 = sprintf('TARGET_STACK="%s";', rc.stack);

str8 = sprintf('/groups/flyTEM/flyTEM/render/bin/manage-stack.sh ${PROJECT_PARAMS} --action CREATE --stack ${TARGET_STACK};');
strcmd = [str1 str3 str8];

try
    [a, resp] = system(strcmd);
catch err_cmd_exec
    kk_disp_err(err_cmd_exec);
    error(['Error executing: ' strcmd]);
end

if strfind(resp, 'caught exception'),
    disp(resp);
    error('create_renderer_stack: server-side error reported');
end


if verbose,
    disp(strcmd);
    disp(a);
    disp(resp);
end

%%
function check_input(rc)
if ~isfield(rc, 'baseURL'), disp_usage; error('baseURL not provided');end
if ~isfield(rc, 'owner'), disp_usage; error('owner not provided');end
if ~isfield(rc, 'project'), disp_usage; error('project not provided');end
if ~isfield(rc, 'stack'), disp_usage; error('stack not provided');end


%%
function disp_usage()
disp('Usage:');
disp('Provide an input struct with fields: baseURL, owner, project, stack');