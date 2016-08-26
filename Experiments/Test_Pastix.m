%% solve Ax=b using Pastix
clc
num_slots = 2;
executable_cmd = sprintf('/groups/flyTEM/home/khairyk/downloads/pastix_gennady/align_tiles_pastix');
lin_system_fn = '/nobackup/flyTEM/khairy/FAFB00v13/matlab_large_matrix_system/two_tile_matrix_system_Ab.mat';
mpi_cmd = sprintf('mpirun -np %d %s %s;', num_slots, executable_cmd, lin_system_fn);
str_env = 'source /groups/flyTEM/home/khairyk/downloads/pastix_gennady/bashrc_for_pastix_20160824';
%str_env = 'export I_MPI_FABRICS=shm:ofa;export I_MPI_FALLBACK=disable;export LD_LIBRARY_PATH=$HOME/lib64:/usr/local/hwloc-1.11.3/lib:$LD_LIBRARY_PATH';





str = sprintf('qsub -pe impi %d "%s;%s"', num_slots, str_env, mpi_cmd);

%str = sprintf('qlogin -pe impi %d -now no;%s', num_slots, mpi_cmd);
%s;qsub -l short=true -N %s %s %s -cwd -V -b y -pe batch 1 "%s" 



%%

% % 
% % 
% % All the you need in order to compile and run the PaStiX-based executable, align_tiles_pastix, is stored in the folder:
% % 
% % 
% %     PASTIX_HOME=/tier2/flyTEM/denisovg/pastix
% % 
% % 
% % Steps to proceed:
% % 
% % 
% % 1) make use of the provided .bashrc file:
% % 
% % 
% %     cp $PASTIX_HOME/bashrc_for_pastix_20160824 ~/.bashrc
% % 
% %     source ~/.bashrc
% % 
% % 
% % 2) qlogin to cluster:
% % 
% % 
% %     qlogin -pe impi <num_slots> -now n
% % 
% % 
% % where <num_slots> is the number of the worker nodes you plan to use (say, = 10)
% % 
% % 
% % 3) compile the executable:
% % 
% % 
% %     cd $PASTIX_HOME
% % 
% %     source compile_align_tiles.sh
% % 
% % 
% % 4) run the executable using the specified <num_slots> on the provided input file:
% % 
% % 
% %     mpirun -np <num_slots> ./align_tiles_pastix ntiles34_axb.mat
% % 
