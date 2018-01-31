function [PM, zdel, c, tId_del] = system_solve_helper_point_match_connectivity(nfirst, nlast, rc, pm, opts)

zdel = [];
tId_del = [];
    if ~isfield(opts, 'filter_point_matches'), opts.filter_point_matches = 1;end
    if ~isfield(opts, 'use_peg'), opts.use_peg = 0;end
    if ~isfield(opts, 'nbrs_step'), opts.nbrs_step = 1;end
    

    dir_scratch = [opts.dir_scratch '/temp_' num2str(randi(3000000))];
    kk_mkdir(dir_scratch);
    % obtain actual section zvalues in given range their ids and also of possible reacquires
    [zu, sID, sectionId, z, ns] = get_section_ids(rc, nfirst, nlast);

    %disp('Loading transformations and tile/canvas ids from Renderer database.....');
    [T, map_id, tIds, z_val] = load_all_transformations(rc, zu, dir_scratch);
    ntiles = size(T,1);

    %% Load point-matches

% disp('loading point matches');

    [M, adj, W, np, discard] = system_solve_helper_load_point_matches(...
        zu, opts, pm, map_id, sID, size(T,1));
    PM.M = M;
    PM.adj = adj;
    PM.W = W;
    PM.np = np;


%%% test for number of connected components (should be one)
%%% could take up a lot of time/memory to do for large PM structs
indx = sub2ind([ntiles ntiles], PM.adj(:,1), PM.adj(:,2));
indx = [indx; sub2ind([ntiles ntiles], PM.adj(:,2), PM.adj(:,1))];
A = sparse(ntiles,ntiles);
A(indx) = 1;
G = graph(A);
c = conncomp(G);
nbins = max(c);
disp([num2str(nfirst) ' to ' num2str(nlast) ': Number of connected components is ' num2str(nbins)]);
if nbins>1
    warning(' Number of connected components in point-match graph cannot exceed one');
    disp('Affected sections:');
    zdel = z_val(find(diff(c)>0));
    disp(zdel);
 tId_del = tIds(c>1);
% %     %%%% query whether user wants to delete those sections or not
% % %     str = input(' --------   Do you want to delete these sections in your source collection? Y/N [N] ', 's');
% % %     if strcmp(str,'Y')
% %         disp(['Deleting ' num2str(numel(zdel)) ' sections.']);
% %         %%%% delete the sections
% %         parfor ix = 1:numel(zdel)
% %             resp = delete_renderer_section(rc, zdel(ix), 0);
% %         end
% %         set_renderer_stack_state_complete(rc);
% %         %%%%
% %         disp('Re-running point-match loading and connectivity analysis');
% %         [PM, zdel, c] = system_solve_helper_point_match_connectivity(nfirst, nlast, rc, pm, opts);
% % %     end
else
    %disp('Found one connected component!');
end
%% 
