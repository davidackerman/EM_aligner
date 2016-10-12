function [mL,err,R, A, b, B, d, Wmx, K, Lm, xout, LL2, U2, tB, td, invalid, PM] = ...
         system_solve(nfirst, nlast, rc, pm, opts, rcout)

%% read tiles, load point-matches, generate linear system and solve
mL = [];
err = [];
R = [];
A = [];
b = [];
B = [];
d = [];
Wmx = [];
K= [];
Lm = [];
xout = [];
LL2 = [];
U2 = [];
tB = [];
td = [];
invalid = [];
PM = [];

dir_scratch = [opts.dir_scratch '/temp_' num2str(randi(3000000))];
kk_mkdir(dir_scratch);
cd(dir_scratch);
diary on;
% obtain actual section zvalues in given range their ids and also of possible reacquires
[zu, sID, sectionId, z, ns] = get_section_ids(rc, nfirst, nlast);

% load all tiles in this range and pool into Msection object
disp('Loading transformations and tile/canvas ids from Renderer database.....');
[T, map_id, tIds, z_val] = load_all_transformations(rc, zu, dir_scratch);
ntiles = size(T,1);
disp(['..system has ' num2str(ntiles) 'tiles...']);
%[L, map_id, tIds] = load_all_tiles(rc,zu);ntiles = numel(L.tiles);
degree = opts.degree;
tdim = (degree + 1) * (degree + 2)/2; % number of coefficients for a particular polynomial
tdim = tdim * 2;        % because we have two dimensions, u and v.
ncoeff = ntiles*tdim;

disp('....done!');diary off;diary on;
%% predict sequence of PM requests to match sequence required for matrix A

sID_all = {};
fac = [];
ismontage = [];
count  = 1;
for ix = 1:numel(zu)   % loop over sections
    %disp(['Montage: ' sID{ix}]);
    sID_all{count,1} = sID{ix};
    sID_all{count,2} = sID{ix};
    ismontage(count) = 1;
    fac(count) = 1;
    count = count + 1;
    for nix = 1:opts.nbrs   % loop over neighboring sections
        if (ix+nix)<=numel(zu)
            %disp(['cross-layer: ' num2str(ix) ' ' sID{ix} ' -- ' num2str(nix) ' ' sID{ix+nix}]);
            sID_all{count,1} = sID{ix};
            sID_all{count,2} = sID{ix+nix};
            ismontage(count) = 0;
            fac(count) = 1/(nix+1);
            count = count + 1;
        end
    end
end
%% perform pm requests
disp('Loading point-matches from point-match database ....');
wopts = weboptions;
wopts.Timeout = 20;
M   = {};
adj = {};
W   = {};
np = {};  % store a vector with number of points in point-matches (so we don't need to loop again later)
parfor ix = 1:size(sID_all,1)   % loop over sections
    %disp([sID_all{ix,1}{1} ' ' sID_all{ix,2}{1} ' ' num2str(ismontage(ix))]);
    if ismontage(ix)
        [m, a, w, n] = load_montage_pm(pm, sID_all{ix,1}, map_id,...
            opts.min_points, opts.max_points, wopts);
    else
        [m, a, w, n] = load_cross_section_pm(pm, sID_all{ix,1}, sID_all{ix,2}, ...
            map_id, opts.min_points, opts.max_points, wopts, fac(ix));
    end
    M(ix) = {m};
    adj(ix) = {a};
    W(ix) = {w};
    np(ix) = {n};  
end
disp('... concatenating point matches ...');
% concatenate
M = vertcat(M{:});
adj = vertcat(adj{:});
W   = vertcat(W{:});
np  = [np{:}]';

% cd(dir_scratch)
% save PM M adj W -v7.3;
% fn = [dir_scratch '/PM.mat'];
% PM = matfile(fn);

disp(' ..... done!');diary off;diary on;
%% generate row slabs of matrix A
disp('Generating system matrix .... ');
split = opts.distribute_A;

npm = size(np,1);
disp(' .... determine row positions of point-pairs (needed for generation of A)...');
n = 2*sum(np);
r_sum_vec = [1;cumsum(2*np(1:npm-1))+1];
pm_per_worker = round(npm/split);
disp([' .... pm_per_worker=' num2str(pm_per_worker)]);
r = zeros(split,2);
for ix=1:split
    pm_min = 1 + (ix-1)*pm_per_worker;
    if ix < split
        pm_max = pm_min   + pm_per_worker-1;
    else
        pm_max = npm;
    end
    r(ix,:) = [pm_min pm_max];
end
indx = find(r(:,1)>npm);
r(indx,:) = [];
r(end,2)  = npm;
split = size(r,1);
disp(' .... export temporary files split_PM_*.mat...');
fn_split = cell(split,1);
for ix = 1:split
    fn_split{ix} = [dir_scratch '/split_PM_' num2str(nfirst) '_' num2str(nlast) '_' num2str(randi(10000000)) '_' num2str(ix) '.mat'];
    vec = r(ix,1):r(ix,2);
    m = M(vec,:);
    a = adj(vec,:);
    ww = W(vec);
    save(fn_split{ix}, 'm', 'a', 'ww');
end
clear M adj W
diary off;diary on;
disp(' .... generate matrix slabs');
degree = 1;
I = {};
J = {};
S = {};
w = {};
parfor ix = 1:split
    [I{ix}, J{ix}, S{ix}, wout] = gen_A_b_row_range(fn_split{ix}, ...
        degree, np,r_sum_vec, r(ix,1), r(ix,2));
    wout(wout==0)= [];
    w{ix} = wout;
end
% delete/cleanup
for ix = 1:split
    try
    delete(fn_split{ix});
    catch err_delete
        kk_disp_err(err_delete);
    end
end


%% collect matrix slabs into one matrix A
disp('.... collect: generate the sparse matrix from I, J and S');
I1 = cell2mat(I(:));clear I;
J1 = cell2mat(J(:));clear J;
S1 = cell2mat(S(:));clear S;
disp('..... done!');
%% solve
disp('Solving ....'); diary off;diary on;
% build system and solve it
A = sparse(I1,J1,S1, n,ntiles*tdim); clear I1 J1 S1;
b = sparse(size(A,1), 1);
w = cell2mat(w(:));
Wmx = spdiags(w,0,size(A,1),size(A,1));
clear w;
d = reshape(T', ncoeff,1);
tB = ones(ncoeff,1);
tB = sparse(1:ncoeff, 1:ncoeff, tB, ncoeff, ncoeff);
K  = A'*Wmx*A + opts.lambda*(tB')*tB;
Lm  = A'*Wmx*b + opts.lambda*(tB')*d;
[x2, R] = solve_AxB(K,Lm, opts, d);
Tout = reshape(x2, tdim, ncoeff/tdim)';% remember, the transformations
disp('.... done!');
%% translate everything to +ve space
delta = 0;
dx = min(Tout(:,3)) + delta;%mL.box(1);
dy = min(Tout(:,6)) + delta;%mL.box(2);

for ix = 1:size(Tout,1)
    Tout(ix,[3 6]) = Tout(ix, [3 6]) - [dx dy];
end
%% ingest into Renderer database
disp('Ingesting data .....');
disp('... export to MET (in preparation to be ingested into the Renderer database)...');

v = 'v1';
if ~stack_exists(rcout)
    disp('.... target collection not found, creating new collection in state: ''Loading''');
    resp = create_renderer_stack(rcout);
end

chks = round(ntiles/32);
cs = 1:chks:ntiles;
cs(end) = ntiles;
disp(' .... ingesting ....');
parfor ix = 1:numel(cs)-1
    vec = cs(ix):cs(ix+1)-1;
    export_to_renderer_database(rcout, rc, dir_scratch, Tout(vec,:),...
                                tIds(vec), z_val(vec), v);
end


% % complete stack
disp(' .... completing stack...');
resp = set_renderer_stack_state_complete(rcout);
disp('.... done!');
diary off;
% delete_renderer_stack(rc_target);






















