classdef Msection
    % Msection  Summary: 
    %           Represents a collection of image tiles that are (partially)
    %           overlapping. These can be one section (all have the same z)
    %           or span mutiple sections, in which case z is that of the
    %           first tile.
    %
    % Msection  Class variables:
    %   z                       - the z coordinate value
    %   tiles                   - array of tile objects
    %   A                       - adjacency matrix (sparse)
    %   method                  - the tile-tile registration method ('alignTEM')
    %
    % Msection  Class Instantiation:
    %           *Scenario 1: obj = Msection();                              %no arguments: empty Msection object
    %           *Scenario 2: obj = Msection(tile_array);                    %tile_array: matlab array of tile objects
    %           *Scenario 3: obj = Msection(layout_file, z);                %full layout file path, and 'z' section
    %
    % Msection  Example usage:
    %           >fn = '/scratch/myexperiments/layoutdata/layout.txt';       % full UNIX path to layout file
    %           >z = 2299;                                                  % select a 'z' section
    %           >L = Msection(fn,z);                                        % create the Msection object with tile information from 'z'
    %           >rL = register(L);                                          % perform registration (stitching) of tiles using currently configured   method
    %           >report(rL);                                                % generate a report summary of registration process
    %           >layer_explorer(rL);                                        % explore registration feedback, errors and perform rudementary editing
    %
    % Msection  Class methods (public):
    %           register
    %           report
    %           show_map                                                    % shows cartoon (rectangles) of the section layout
    %           show_map_tile                                               % shows neighborhood (connections) of a particular tile
    %           show_tile                                                   % shows the image of the tile and its neighbors overlayed
    %           mi_calc                                                     % calculates the mutual information criterion for each adjacent pair of tiles
    %           add_tile
    %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        z;                              % the z id (numeric a la Karsh) as it appears in the layout file
        stackID = 0;                    % identifying the section as belonging to stack with same ID
        sectionID;                      % Unique ID appearing in the renderer databasse
        tiles;                          % array of tiles
        G;
        A;                              % adjacency matrix
        original_layout_file = [];      % fully qualified unix path to layout file
        method = 'alignTEM';            % use alignTEM alignment method as default
        viewer = 'matlab';            % select which viewer to use to look at layer -- 'renderer' or 'matlab' ;
        regConf = [];                   % struct for the registration method configuration
        rprt = [];                      % struct to store all registration reports
        tile_show_limit = 20;           % limit of number of tiles to open when show_tile is called
        tile_show_method = 'overlay';   % method used to display the tiles -- default: overlay
        update_adjacency_switch = 1;    % switch indicating 1--> update adjacency matrix, 0 --> not needed
        update_tile_info_switch = -1;   % 1: update each tile individually one time, 0: no need to update, -1: update all to the data in the first tile
        map_display_fac = 1;            % size factor for rectangles representing tiles
        map_id;                         % maps: 'id' to index into the 'tiles' array
        map_renderer_id;                % maps: 'renderer_id' to index into 'tiles' array
        mosaic = [];                    % mosaic image of all tiles
        T = [1 0 0;0 1 0;0 0 1];        % generic identity transform 
        X, Y, box;                      % housekeeping arrays
        pm;                             % point-matches struct
        n_points_pair_same = 10;        % number of points to generate for each image pair in a forward transformation (used by get_point_pairs);
        dthresh_factor = 0.9;           % factor for adjacency calculation. Adjacent tiles have a center that is diagonal*dthresh_fac removed from the center of the tile
        edge_tiles = [];                % logical vector of length numel(obj.tiles) to indicate edge tiles = 1 or non-edge = 0;
    end

    methods
        
        function obj = Msection(arg1, arg2, arg3)   % constructor
            % *Scenario 1: obj = Msection(); %no arguments: empty Msection object
            % *Scenario 2: obj = Msection(rc, z); %  rc is a struct that defines a renderer collection. z is z-value in the renderer collection
            % *Scenario 3: obj = Msection(tile_array); %  tile_array = array of tile objects
            % *Scenario 4: obj = Msection(layout_original, z); %arg1 = full layout file path, arg2 = z section desired
            
            if nargin~=0
                
                %%%% generate Msection object from array of tile objects
                if nargin==1
                    %%% in this case we are given an array of tiles
                    obj.z = arg1(1).z; % use z value provided in first tile 
                    obj.tiles = arg1;
                    obj.original_layout_file = [];
                    obj = update_adjacency(obj);
                    obj = generate_hash_tables(obj);
                end
                
                %%%% import from renderer database
                if nargin==2 && isstruct(arg1)    % then we construct a section based on renderer information
                    % % Example rc (arg1)
                    % % rc.stack = 'v9_acquire_LC_merged_2';
                    % % rc.owner='flyTEM';
                    % % rc.project='FAFB00';
                    % % rc.server='http://tem-services.int.janelia.org:8080/render-ws/v1';
                    rc = arg1;
                    z  = arg2; % is the z position of the section
                    
                    % % get a list of all tiles 
                    urlChar = sprintf('%s/owner/%s/project/%s/stack/%s/z/%d/tile-specs', ...
                        rc.server, rc.owner, rc.project, rc.stack, z);
                    j = webread(urlChar);
                    % generate the tiles
                    jt = tile;
                    sectionID = j(1).layout.sectionId;
                    disp('Loading section ... ');
                    parfor jix = 1:numel(j)
                        jt(jix) = tile(j(jix));
                        jt(jix).z = z;
                    end
                    % loop over tiles to set tile id
                    parfor ix = 1:numel(jt)
                        jt(ix).id = ix;
                    end
                    obj = Msection(jt);
                    obj = update_tile_sources(obj, rc.owner, rc.project, rc.stack, rc.server);
                    obj.sectionID = sectionID;
                end
                
                %%%%%% import from layout file
                if nargin==3 && iscell(arg1) && strcmp(arg2, 'BK_format')
                    %% A layer is constructed when called from the Mstack generate_stack function.
                    %%% arg1 is a cell array containing all information for
                    %%% that layer, and arg2 indicates it is in BK_format
                    %%% (i.e. the data was read from BK format layout file.
                    %%% This also means we are using affine2d
                    %%% transformations
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    Xo = 0; % translate all tiles by that X and Y while ingesting
                    Yo = 0;
                    id_vec = cell(numel(arg1{1}), 1);
                    t = tile;
                    X = 0;
                    Y = 0;
                    for tix = 1:numel(arg1{1})      %% we have to generate the individual tile objects, so let's loop over them
                       
                        if numel(arg1)==12
                        t(tix) = tile(arg1{1}(tix), arg1{2}(tix), arg1{3}(tix), arg1{4}(tix), ...
                            arg1{5}(tix)-Xo, arg1{6}(tix), arg1{7}(tix), arg1{8}(tix)-Yo, arg1{9}(tix), ...
                            arg1{10}(tix), arg1{11}(tix), arg1{12}{tix});
                        elseif numel(arg1)==13
                            t(tix) = tile(arg1{1}(tix), arg1{2}(tix), arg1{3}(tix), arg1{4}(tix), ...
                            arg1{5}(tix)-Xo, arg1{6}(tix), arg1{7}(tix), arg1{8}(tix)-Yo, arg1{9}(tix), ...
                            arg1{10}(tix), arg1{11}(tix), arg1{12}{tix}, arg1{13}(tix));
                        elseif numel(arg1)==14
                            t(tix) = tile(arg1{1}(tix), arg1{2}(tix), arg1{3}(tix), arg1{4}(tix), ...
                                arg1{5}(tix)-Xo, arg1{6}(tix), arg1{7}(tix), arg1{8}(tix)-Yo, arg1{9}(tix), ...
                                arg1{10}(tix), arg1{11}(tix), arg1{12}{tix}, arg1{13}(tix), arg1{14}(tix));
                        elseif numel(arg1)==15
                            t(tix) = tile(arg1{1}(tix), arg1{2}(tix), arg1{3}(tix), arg1{4}(tix), ...
                                arg1{5}(tix)-Xo, arg1{6}(tix), arg1{7}(tix), arg1{8}(tix)-Yo, arg1{9}(tix), ...
                                arg1{10}(tix), arg1{11}(tix), arg1{12}{tix}, arg1{13}(tix), arg1{14}(tix), arg1{15}{tix});
                        end
                        
                        id_vec{tix} = arg1{2}(tix);
                        X(tix) = arg1{5}(tix);
                        Y(tix) = arg1{8}(tix);
                    end
                    obj.tiles = t;
                    obj.X = X;
                    obj.Y = Y;
                    obj.z = arg1{1}(1);
                    obj.original_layout_file = [];
                    obj.map_id = containers.Map(id_vec, 1:numel(arg1{1}));     % generate a hash table with id's as the keys
                    obj = update_adjacency(obj);
                    obj = generate_hash_tables(obj);
                end
                
                %%%%%% import from layout file
                if nargin==2 && ischar(arg1)
                    %% in this case we construct a layer by directly reading the relevant z tile information from the layout file
                    obj = Msection;
                    obj = import_from_layout_txt(obj, arg1, arg2);
                    obj = generate_hash_tables(obj);
                    obj = update_adjacency(obj);
                end
                
                %%%%%% copy -- generate a section based on another
                if nargin ==3 && strcmp('Msection', class(arg1))
                    % then we are generating a subsection from a section
                    % arg1 is the section, arg2 the col_range and arg3 the
                    % row_range
                    obj = Msection;
                    obj = get_subsection(obj, arg1, arg2, arg3);
                    obj = generate_hash_tables(obj);
                    obj = update_adjacency(obj);
                end
            end
            if ~isempty(obj.tiles), obj = update_tile_info(obj);end
            obj = configure_registration(obj);
            obj = tile_pos_gen(obj);
            %if nnz(obj.A)==0, warning('No adjacency information in A');end
        end
    end
end
























