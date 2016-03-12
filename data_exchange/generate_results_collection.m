function generate_results_collection(rc)
% INPUT: 
% rc is a struct that provides renderer database configuration
% rc = 
% 
%       owner: 'khaled'
%     project: 'test01'
%       source_stack: 'stack01'
%   target_stack:
%   CATMAID_ROOT_DIR:
%     baseURL: 'http://tem-services:8080/render-ws/v1'
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(rc, 'verbose')
    verbose = rc.verbose;
else
    verbose = 1;
end
% %%%%%%%%%%%%%%%%%%%  remove collection of same name if already exists
str1 = sprintf('PROJECT_PARAMS="--baseDataUrl %s --owner %s --project %s";', rc.baseURL, rc.owner, rc.target_project);
str2 = sprintf('TARGET_STACK="%s";', rc.target_stack);
str3 = sprintf('/groups/flyTEM/flyTEM/render/bin/manage-stack.sh ${PROJECT_PARAMS} --action DELETE --stack ${TARGET_STACK}');
strcmd = [str1 str2 str3];
[a, resp] = system(strcmd);
if verbose,
disp(a);
disp(resp);
end
% #################### create
if verbose, disp('------------------------ Configure --------------');end

str1_source = sprintf('PROJECT_PARAMS="--baseDataUrl %s --owner %s --project %s";', rc.baseURL, rc.owner, rc.source_project);  
str1_target = sprintf('PROJECT_PARAMS="--baseDataUrl %s --owner %s --project %s";', rc.baseURL, rc.owner, rc.target_project); 
target_project = sprintf('TARGET_PROJECT="%s";', rc.target_project);
str2 = sprintf('SOURCE_STACK="%s";', rc.source_stack);
str3 = sprintf('TARGET_STACK="%s";', rc.target_stack);
str4 = sprintf('OWNER="%s";', rc.owner);
str5 = sprintf('PROJECT="%s";', rc.project);
str6 = sprintf('MET_FORMAT="%s";', rc.MET_format);
%str7 = sprintf('CATMAID_ROOT_DIR="%s";', rc.CATMAID_ROOT_DIR);



% # MET Formats:
% #   "v1" is the original format that expects all tiles to have 6 affine parameters
% #   "v3" is the new polynomial format where the number of parameters is identified for each tile 
% MET_FORMAT="v3"
% CATMAID_ROOT_DIR="/tier2/flyTEM/nobackup/rendered_boxes/low_latency_montage"

%%%%%% create stack
if verbose, disp(' ----------------Creating stack---------------');end
str8 = sprintf('/groups/flyTEM/flyTEM/render/bin/manage-stack.sh ${PROJECT_PARAMS} --action CREATE --stack ${TARGET_STACK};');
strcmd = [str1_target str3 str8];
[a, resp] = system(strcmd);
if verbose,
disp(strcmd);
disp(a);
disp(resp);
end
if strfind(resp, 'caught exception'), error('import error detected');end



% load MET files
if verbose, disp(' ----------- Loading MET files --------');end
str9 = sprintf('MEMORY="1G";');
str10 = sprintf('JAVA_CLASS="org.janelia.render.client.ImportMETClient";');

% <sosi --------------------- maybe loop over list of MET files when doing
% slab in the future  --------------- /sosi>
MET_file = rc.MET_files{1};
str11 = sprintf('MET_FILE="%s";', MET_file);
str12 = sprintf('/groups/flyTEM/flyTEM/render/pipeline/bin/run_ws_client.sh ${MEMORY} ${JAVA_CLASS} ${PROJECT_PARAMS} --targetProject ${TARGET_PROJECT} --stack ${SOURCE_STACK} --targetStack ${TARGET_STACK} --metFile ${MET_FILE} --formatVersion ${MET_FORMAT};');
strcmd = [str1_source target_project str2 str3 str6 str9 str10 str11 str12];
[a, resp] = system(strcmd);
if verbose
disp(strcmd);
disp(a);
disp(resp);
end
if strfind(resp, 'caught exception'), error('import error detected');end


% # complete montage stack
if verbose, disp(' ---------- Completing stack -----------');end
str13 = sprintf('/groups/flyTEM/flyTEM/render/bin/manage-stack.sh ${PROJECT_PARAMS} --action SET_STATE --stackState COMPLETE --stack ${TARGET_STACK}');
strcmd = [str1_target str3 str13];
[a, resp] = system(strcmd);
if verbose,
disp(strcmd);
disp(a);
disp(resp);
end
if strfind(resp, 'caught exception'), error('import error detected');end










%% if isfield(rc, 'verbose')
%     verbose = rc.verbose;
% else
%     verbose = 1;
% end
% % %%%%%%%%%%%%%%%%%%%  remove collection of same name if already exists
% str1 = sprintf('PROJECT_PARAMS="--baseDataUrl %s --owner %s --project %s";', rc.baseURL, rc.owner, rc.project);
% str2 = sprintf('TARGET_STACK="%s";', rc.target_stack);
% str3 = sprintf('/groups/flyTEM/flyTEM/render/bin/manage-stack.sh ${PROJECT_PARAMS} --action DELETE --stack ${TARGET_STACK}');
% strcmd = [str1 str2 str3];
% [a, resp] = system(strcmd);
% if verbose,
% disp(a);
% disp(resp);
% end
% % #################### create
% if verbose, disp('------------------------ Configure --------------');end
% 
% str1 = sprintf('PROJECT_PARAMS="--baseDataUrl %s --owner %s --project %s";', rc.baseURL, rc.owner, rc.project);
% str2 = sprintf('SOURCE_STACK="%s";', rc.source_stack);
% str3 = sprintf('TARGET_STACK="%s";', rc.target_stack);
% str4 = sprintf('OWNER="%s";', rc.owner);
% str5 = sprintf('PROJECT="%s";', rc.project);
% str6 = sprintf('MET_FORMAT="%s";', rc.MET_format);
% %str7 = sprintf('CATMAID_ROOT_DIR="%s";', rc.CATMAID_ROOT_DIR);
% 
% 
% 
% % # MET Formats:
% % #   "v1" is the original format that expects all tiles to have 6 affine parameters
% % #   "v3" is the new polynomial format where the number of parameters is identified for each tile 
% % MET_FORMAT="v3"
% % CATMAID_ROOT_DIR="/tier2/flyTEM/nobackup/rendered_boxes/low_latency_montage"
% 
% %%%%%% create montage stack
% if verbose, disp(' ----------------Creating stack---------------');end
% str8 = sprintf('/groups/flyTEM/flyTEM/render/bin/manage-stack.sh ${PROJECT_PARAMS} --action CREATE --stack ${TARGET_STACK};');
% strcmd = [str1 str3 str8];
% [a, resp] = system(strcmd);
% if verbose,
% disp(strcmd);
% disp(a);
% disp(resp);
% end
% if strfind(resp, 'caught exception'), error('import error detected');end
% 
% 
% 
% % load MET files
% if verbose, disp(' ----------- Loading MET files --------');end
% str9 = sprintf('MEMORY="1G";');
% str10 = sprintf('JAVA_CLASS="org.janelia.render.client.ImportMETClient";');
% 
% % <sosi --------------------- maybe loop over list of MET files when doing
% % slab in the future  --------------- /sosi>
% MET_file = rc.MET_files{1};
% str11 = sprintf('MET_FILE="%s";', MET_file);
% str12 = sprintf('/groups/flyTEM/flyTEM/render/pipeline/bin/run_ws_client.sh ${MEMORY} ${JAVA_CLASS} ${PROJECT_PARAMS} --acquireStack ${SOURCE_STACK} --alignStack ${TARGET_STACK} --metFile ${MET_FILE} --formatVersion ${MET_FORMAT};');
% strcmd = [str1 str2 str3 str6 str9 str10 str11 str12];
% [a, resp] = system(strcmd);
% if verbose
% disp(strcmd);
% disp(a);
% disp(resp);
% end
% if strfind(resp, 'caught exception'), error('import error detected');end
% 
% 
% % # complete montage stack
% if verbose, disp(' ---------- Completing stack -----------');end
% str13 = sprintf('/groups/flyTEM/flyTEM/render/bin/manage-stack.sh ${PROJECT_PARAMS} --action SET_STATE --stackState COMPLETE --stack ${TARGET_STACK}');
% strcmd = [str1 str3 str13];
% [a, resp] = system(strcmd);
% if verbose,
% disp(strcmd);
% disp(a);
% disp(resp);
% end
% if strfind(resp, 'caught exception'), error('import error detected');end







































