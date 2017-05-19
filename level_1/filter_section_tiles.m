function filter_section_tiles(rc, rcsource, z, dir_temp, lambda)



if nargin<4
    lambda = 1.0;
end

verbose = 1;
translate_to_positive_space = 0;
complete = 0;
disableValidation = 0;

for zix = 1:numel(z)
    L = Msection(rc, z(zix)); % read the section with z-value z
    [L] = filter_based_on_tile_area(L, lambda);
    indx = [L.tiles(:).state];
    L.tiles(find(indx==-3))=[];
    delete_renderer_section(rc, z(zix), 0);
    ingest_section_into_renderer_database(L, rc, rcsource, dir_temp, ...
        translate_to_positive_space, complete, disableValidation);
end