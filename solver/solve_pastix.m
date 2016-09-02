function [x2, R] = solve_pastix(A, b,ncpus)
%%% use pastix to solve this system
kk_clock;
tic;
%% configure
disp('Solving Ax=b system using pastix...');
dir_temp_mx = '/nobackup/flyTEM/khairy/FAFB00v13/matlab_large_matrix_system';
%PASTIX_HOME = '/nobackup/flyTEM/khairy/FAFB00v13/pastix_solver';
PASTIX_HOME = '/tier2/flyTEM/denisovg/pastix';

PASTIX_DATA = [dir_temp_mx];
lin_system_fn = [PASTIX_DATA '/linear_system_Ab.mat'];
str_pastix_solver = PASTIX_HOME;
str_pastix_run = [PASTIX_HOME '/run_align_tiles_mod01.sh'];
pastix_script = [str_pastix_solver '/align_tiles_mod01.sh'];
str_env = sprintf('export PASTIX_HOME=%s;export PASTIX_DATA=%s;',str_pastix_solver, PASTIX_DATA); 

kk_mkdir(PASTIX_DATA);



%% save A and b to disk
save(lin_system_fn, 'A', 'b');

%% system command to solve
str = sprintf('%s %d %s %s', pastix_script, ncpus, lin_system_fn, str_pastix_run); 
disp(str);
[a, resp_str] = system(str);disp(resp_str);
[a, resp_str] = system('qstat');disp(resp_str);
%% wait for response
disp('... waiting for solution ...');
dt = 10;
pause(dt);
solving = 1;
while(solving)
    fn = dir(PASTIX_DATA);
    for fix = 2:numel(fn)
        if strcmp(fn(fix).name(1), 'x')
            solving = 0;
        end
    end
    pause(dt);
end
disp('Done.');
%% read and assemble result into x2
%x2 = assemble_x(PASTIX_DATA);
disp('Assembling solution vector...');
fn = dir(PASTIX_DATA);
fn_range1 = [];
fn_range2 = [];%zeros(numel(fn)-2,1);
count = 1;
for fix = 2:numel(fn)
    if strcmp(fn(fix).name(1), 'x')
        c = strsplit(fn(fix).name(3:end), '-');
        fn_range1(count) = str2double(c{1});
        c = strsplit(c{end}, '.');
        fn_range2(count) = str2double(c{1});
        count = count + 1;
    end
end
[~, ia] = sort(fn_range2);
fn_range1 = fn_range1(ia);
fn_range2 = fn_range2(ia);
% now reassemble everything
x2 = [];
for fix = 1:numel(fn_range1)
    fn = sprintf('%s/x_%s-%s.mat', PASTIX_DATA, num2str(fn_range1(fix)), ...
        num2str(fn_range2(fix)));
    disp(fn);
    c = load(fn);
    x2 = [x2;c.x(:)];
end
disp('Done.');
%% calculate residual
R = A*x2-b;
toc

kk_clock
