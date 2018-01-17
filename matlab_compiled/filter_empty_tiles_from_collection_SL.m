function filter_empty_tiles_from_collection_SL(fn)
% Intended for deployment: solve matrix system using translation based on json input provided by fn


% read json input
sl = loadjson(fileread(fn));

if sl.verbose
    kk_clock();
    disp(['Using input file: ' fn]);
    if numel(sl.rc)==1
        disp(['Collection to be filtered: ' sl.rc]);
    else
        disp(['Input collection: ' sl.rc(1)]);
        disp(['Output collection: ' sl.rc(2)]);
    end
    if isfield(sl, 'zs')
        disp(['z range: ' num2str(min(sl.zs)) '-' num2str(max(sl.zs))]);
    end
    if isfield(sl, 'opts')
        disp(['Options: ' sl.opts]);
    end
end

fft2_results = filter_empty_tiles_from_collection(sl.rcs, sl.zs, sl.opts);