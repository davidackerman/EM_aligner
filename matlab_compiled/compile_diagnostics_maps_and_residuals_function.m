%% compile (you must customize this script to your system/environment


cd /groups/flyTEM/home/ackermand/EM_aligner/matlab_compiled

astr = [];

fn = dir('/groups/flyTEM/home/ackermand/EM_aligner/external/jsonlab/*.m');
for ix = 1:numel(fn)
astr = [astr sprintf(' -a /groups/flyTEM/home/ackermand/EM_aligner/external/jsonlab/%s',fn(ix).name)];
end


str = sprintf('mcc -m -R -nodesktop -v diagnostics_maps_and_residuals_function.m %s;', astr);
eval(str);
