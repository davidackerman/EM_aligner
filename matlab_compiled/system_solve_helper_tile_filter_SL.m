function system_solve_helper_tile_filter_SL(fn)
% Intended for deployment: solve matrix system using translation based on json input provided by fn


% read json input
sl = loadjson(fileread(fn));

if sl.verbose
    kk_clock();
    disp(['Using input file: ' fn]);
    disp(['Input collection: ' sl.rcin]);
    disp(['Output collection: ' sl.rcout]);
    if isfield(sl, 'zs')
        if esmpty(sl.zs)
            disp("Will filter over collection's enitre z-range");
        else
            disp(['z range: ' num2str(min(sl.zs)) '-' num2str(max(sl.zs))]);
        end
    end
    disp(['Options: ' sl.opts]);
end
num_cores = feature('numcores');
parpool(num_cores);
current_pool = gcp;
if num_cores ~= current_pool.NumWorkers
    error('The full number of possible cores (%d) are not being used; only %d are being used', num_cores, current_pool.NumWorkers);
end
fft2_results = filter_empty_tiles_from_collection( sl.rcin, sl.rcout, sl.zs, sl.opts, sl.pm );