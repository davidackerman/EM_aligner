function [ output_args ] = remove_blank_tiles_from_collection(rc, zu, opts )
if nargin==1
    zu = get_section_ids(rc);
    opts = [];
end
if nargin==2
    if isstruct(varargin{1})
        opts = varargin{1};
        zu = get_section_ids(rc);
    else
        opts = [];
        zu = varargin{1};
    end
end
if nargin==3
    zu = varargin{1};
    opts = varargin{2};
end

if ~isfield(opts, 'cutoff_value'), opts.cutoff_value = 1E4; end
if ~isfield(opts, 'save_removed_tile_images'), opts.save_removed_tile_images = 1;  end
if ~isfield(opts, 'output_directory'), opts.output_directory = pwd;
fft2_results  = find_blank_tiles( rc, zu, opts );
all_tile_ids = [fft2_results(:).tile_ids];
all_tile_metrics = [fft2_results(:).blank_metric];
all_tile_is_in_middle = [fft2_results(:).is_in_middle];
all_tile_zs = [fft2_results(:).tile_zs];
tile_indices_below_cutoff = all_tile_metrics<opts.cutoff_value;
tile_ids_below_cutoff= all_tile_ids(tile_indices_below_cutoff);
tile_indices_to_replace = find(all_tile_is_in_middle & tile_indices_below_cutoff);
tile_indices_to_remove = find(~all_tile_is_in_middle & tile_indices_below_cutoff);
fid_remove = fopen([opts.output_directory '/removed_tiles.txt'],'wt');
for i=1:numel(tile_indices_to_remove)
   fprintf(fid_remove, '%s \n', all_tile_ids{tile_indices_to_remove(i)});
end
fclose(fid_remove);
% replace middle tiles with neighboring ones
fid_replace = fopen([opts.output_directory '/replaced_tiles.txt'],'wt');
for i=1:numel(tile_indices_to_replace)
    %try previous section first
    bad_tile_z = all_tile_zs(tile_indices_to_replace(i));
    L_bad = Msection(rc, bad_tile_z);
    bad_tile_index = find(ismember({L_bad.tiles.renderer_id}, all_tile_ids(tile_indices_to_replace(i))));
    bad_tile_column = L_bad.tiles(bad_tile_index).col;
    found_a_good_tile = 0;
    count = 1;
    while ~found_a_good_tile
        good_tile_z = bad_tile_z - (-1)^count * ceil(count/2); %will search +/- 1 section, 2 sections, 3...
        if good_tile_z>=min(all_tile_zs) && good_tile_z <=max(all_tile_zs)
            L_good = Msection(rc, good_tile_z);
            good_tile_index = find(ismember([L_good.tiles.col],bad_tile_column));
            if ~isempty(good_tile_index) && sum(ismember(tile_ids_below_cutoff, L_good.tiles(good_tile_index).renderer_id))==0 %then the good tile is indeed a good tile
                found_a_good_tile = 1;
            end
        end
        count = count+1;
        if count>10
           error('Too many adjacent blank sections'); 
        end
    end
    path = L_bad.tiles(bad_tile_index).path;
    [pathstr, name, ext] = fileparts(path);
    %system(['mkdir -p ' pathstr '/badImages/']);
    %system(['mv ' pathstr '/' name ext ' ' pathstr '/badImages/ ']);
    %system(['ln -s ' L_good.tiles(good_tile_index).path ' ' L_bad.tiles(bad_tile_index).path]);
    fprintf(fid_replace, '%s %s \n', all_tile_ids{tile_indices_to_replace(i)}, L_good.tiles(good_tile_index).path);
end
fclose(fid_replace);
% remove edge blank tiles
tile_ids_to_remove = all_tile_ids(tile_indices_to_remove);
delete_renderer_tile(tile_ids_to_remove);
end

