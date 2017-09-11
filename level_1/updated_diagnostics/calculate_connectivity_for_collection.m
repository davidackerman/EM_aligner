function [ S, G, tIds, clusters, small_clusters] = calculate_connectivity_for_collection( zs, rc, point_matches, options, min_cluster_size)
if isempty(zs)
    [zu, sID, ~, ~, ~, ~] = get_section_ids(rc);
else
    [zu, sID, ~, ~, ~, ~] = get_section_ids(rc,zs(1),zs(end));
end
[T, map_id, tIds, ~, r, c] = load_all_transformations(rc, zu, options.dir_scratch);
[M, adj, ~, np] = system_solve_helper_load_point_matches(zu, options,point_matches, map_id, sID, size(T,1), r, c);
number_of_pairs = size(M,1);
tiles_1 = zeros(1, number_of_pairs);
tiles_2 = zeros(1, number_of_pairs);
number_of_matches = zeros(1, number_of_pairs);
parfor i=1:number_of_pairs
    tiles_1(i) = adj(i,1);
    tiles_2(i) = adj(i,2);
    number_of_matches(i) = np(i);
end
number_of_tiles = numel(tIds);
S = sparse(tiles_1, tiles_2, number_of_matches, number_of_tiles, number_of_tiles);
G = graph(S,'upper');
clusters = conncomp(G);
[counts,cluster_number] = histcounts(clusters,unique(clusters));
small_clusters = cluster_number(counts<min_cluster_size);
%for i=1:numel(small_clusters)
%    fprintf('%s\n',tIds{clusters==small_clusters(i)});
%end
end

