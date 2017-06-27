function [D, I, x1, y1] = detect_jumps_drift_warps(...
    rc, rc_vec, nfirst, nlast, nprobes, delta, delta_z, scale, x1, y1)
% evaluate alignment by detecting drift, jumps or warps at a specified number
% nspots of spots based on corr2 (in jump_detection.m). this is done for each
% slab of thickness delta_z sections
% D is a cell array of dimensions nslabs x 1  
% nslabs is determined as nfirst:delta_z:nlast
% Also:
% D{1} is of dimensions  (collection number (including reference)) x nprobes
% each element of D{1} is of dimension delta_z x 1 and includes cross 
% correlations within that probe
% I is a cell array (of dimension nslabs) of cell arrays of probes,
% each is an image volume for that probe.



% check that all collections actually exist
res = stack_exists(rc); if res==0, disp(rc);error('Reference stack does not exist');end
for rcix = 1:numel(rc_vec)
    res = stack_exists(rc_vec(rcix)); 
    if res==0, disp(rc_vec(rcix));
        error('Test stack does not exist');end
end


%% determine slabs
slabs = nfirst:delta_z:nlast;
disp(slabs);
D = cell(numel(slabs)-1,1);
I = cell(numel(slabs)-1,numel(rc_vec)+1);
for six = 1:numel(slabs)-1
    try
    disp(['Processing slab ' num2str(six) ' of ' num2str(numel(slabs)-1)]);
    %  provide a set of centers for boxes for current slab
    z_first = slabs(six);
    z_last = slabs(six+1);
    if nargin<9
    disp('---Obtaining candidate probe centers');
    [x1, y1] = get_tile_centers(rc, z_first , 0);  % just use the set of tile centers
    indx = randi(size(x1,1), nprobes,1); % determine indices of centers to use
    x1 = x1(indx);
    y1 = y1(indx);
    else
        x1 = x1(:);
        y1 = y1(:);
        nprobes = size(x1,1);
        if ~(size(x1,1)==size(y1,1)), error('Input x and y must be equal in size');end
    end
    
    disp('-----Processing reference collection:'); 
    disp(rc);
    [d(1,:), im] = jump_detection(rc, z_first, z_last, delta, scale, [x1 y1]);
    I{six,1} = im; %% cell array of image volumes corresponding to number of probes (size(x1,1))
    % loop over test Renderer collections
    for rcix = 1:numel(rc_vec)
        disp('-----Processing test collection:');
        disp(rc_vec(rcix));
        % convert reference [x1 y1] (in world coordinates of reference stack rc)
        % to world coordinates of the next
        x = zeros(nprobes,1);
        y = zeros(nprobes,1);
        for pix = 1:nprobes
            try
            [x(pix), y(pix)] = world_to_world(rc, rc_vec(rcix), x1(pix), y1(pix), z_first);
            catch err_point_not_found
               x(pix) = nan;
               y(pix) = nan;
            end
        end
        
        [d(rcix+1,:), im] = jump_detection(rc_vec(rcix), z_first, z_last, delta, scale, [x y]);
        I{six,rcix+1} = im;
    end
    D{six} = d;
    catch err_slab_processing
        disp('Error ---- caught exception:');
        disp(six);
       disp_kk_err(err_slab_processing);
    end
    
end
