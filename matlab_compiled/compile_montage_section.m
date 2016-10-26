%% compile (you must customize this script to your system/environment

dir_curr = pwd;

cd /groups/flyTEM/home/khairyk/EM_aligner/matlab_compiled

astr = [];
fn_use = dir('/groups/flyTEM/home/khairyk/EM_aligner/classes/@Msection*.m');
for ix = 1:numel(fn_use)
astr = [astr sprintf(' -a /groups/flyTEM/home/khairyk/EM_aligner/classes/@Msection/%s',fn_use(ix).name)];
end

%astr = [];
fn_use = dir('/groups/flyTEM/home/khairyk/EM_aligner/classes/@tile*.m');
for ix = 1:numel(fn_use)
astr = [astr sprintf(' -a /groups/flyTEM/home/khairyk/EM_aligner/classes/@tile/%s',fn_use(ix).name)];
end

fn_use = dir('/groups/flyTEM/home/khairyk/EM_aligner/solver/*.m');
for ix = 1:numel(fn_use)
astr = [astr sprintf(' -a /groups/flyTEM/home/khairyk/EM_aligner/solver/%s',fn_use(ix).name)];
end

%astr = [];
fn_use = dir('/groups/flyTEM/home/khairyk/mwork_nogit/fileexchange/jsonlab/*.m');
for ix = 1:numel(fn_use)
astr = [astr sprintf(' -a /groups/flyTEM/home/khairyk/mwork_nogit/fileexchange/jsonlab/%s',fn_use(ix).name)];
end

str_compile = sprintf('mcc -m -R -nodesktop -v montage_section_SL_prll.m %s;', astr);
eval(str_compile);



% cd /groups/flyTEM/home/khairyk/EM_aligner/matlab_compiled
%%
disp('Copy process started ....');

copyfile('/groups/flyTEM/home/khairyk/EM_aligner/matlab_compiled/montage_section_SL_prll', ...
         '/groups/flyTEM/flyTEM/matlab_compiled/sl7/montage_section_SL_prll');
 
copyfile('/groups/flyTEM/home/khairyk/EM_aligner/matlab_compiled/montage_section_SL_prll', ...
         '/gpfs1/scratch/spc/matlab_work/montage/montage_section_SL_prll');
     
copyfile('/groups/flyTEM/home/khairyk/EM_aligner/matlab_compiled/sample_montage_input.json', ...
         '/groups/flyTEM/flyTEM/matlab_compiled/sl7/sample_montage_input.json');

disp('Done!');

cd(dir_curr);