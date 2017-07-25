function [M, adj, W, np, discard] = system_solve_helper_load_point_matches(...
    zu, opts, pm, map_id, sID, ntiles, r, c)
% load all available point matches
% ilter them,
% and after that select points randomly to limit the size of
% point-matches block
if nargin<7  % then we don't use row/col information
    r = [];
    c = [];
end
spmd
    warning('off', 'vision:obsolete:obsoleteFunctionality');
end

disp('** STEP 2:  Load point-matches ....');
disp(' ... predict sequence of PM requests to match sequence required for matrix A');
sID_all = {};
fac = [];
ismontage = [];
count  = 1;
for ix = 1:numel(zu)   % loop over sections  -- can this be made parfor?
    %disp(['Setting up section: ' sID{ix}]);
    sID_all{count,1} = sID{ix};
    sID_all{count,2} = sID{ix};
    ismontage(count) = 1;
    fac(count) = 1;
    count = count + 1;
    for nix = 1:opts.nbrs
        if (ix+nix)<=numel(zu)
            %disp(['cross-layer: ' num2str(ix) ' ' sID{ix} ' -- ' num2str(nix) ' ' sID{ix+nix}]);
            sID_all{count,1} = sID{ix};
            sID_all{count,2} = sID{ix+nix};
            ismontage(count) = 0;
            fac(count) = opts.xs_weight/(nix+1);
            count = count + 1;
        end
    end
end
% clear sID
% % perform pm requests
disp('Loading point-matches from point-match database ....');
wopts = weboptions;
wopts.Timeout = 20;

M   = {};
adj = {};
W   = {};
np = {};  % store a vector with number of points in point-matches (so we don't need to loop again later)
for ix = 1:size(sID_all,1)   % loop over sections
    %disp([sID_all{ix,1}{1} ' ' sID_all{ix,2}{1} ' ' num2str(ismontage(ix))]);
    % when loading point matches, load all available point
    % then filter them, and after that select points randomly to limit the size of
    % point-matches block
    if ismontage(ix)
        [m, a, w, n, xnp] = load_montage_pm(pm, sID_all{ix,1}, map_id,...
            opts.min_points, inf, wopts, r, c);
    else
        [m, a, w, n, xnp] = load_cross_section_pm(pm, sID_all{ix,1}, sID_all{ix,2}, ...
            map_id, opts.min_points, inf, wopts, fac(ix));
    end
    
    
    %%% filter set of point-matches for this pair
    if opts.filter_point_matches
        if sum(n)
            warning off;
            [m, w, a, n] = filter_pm_local(m, w, a, n, opts.pmopts);
            warning on;
        end
    end
    % loop over point-match sets: reduce to max_points if necessary
    delix = [];
    for pmix = 1:size(m,1)    % loop over point-match sets
        pmm = m(pmix,:);
        pmm1 = pmm{1};
        if size(pmm1,1)>opts.max_points  % do we have more than opts.max_points point-matches in the set pmix
            indx = randi(size(pmm1,1)-1, opts.max_points,1);  % define random indices
            pmm2 = pmm{2};
            pmm1 = pmm1(indx,:);
            pmm2 = pmm2(indx,:);
            m(pmix,:) = {pmm1,pmm2};
            pmw = w{pmix};
            w{pmix} = pmw(indx);
            n(pmix) = length(w{pmix});
        end
        
        if size(pmm1,1)<opts.min_points % if we have insufficient points
            delix = [delix pmix];
        end
    end
    if ~isempty(delix)
        m(delix,:) = [];
        w(delix) = [];
        n(delix) = [];
        a(delix,:) = [];
    end
    
    M(ix) = {m};
    adj(ix) = {a};
    W(ix) = {w};
    np(ix) = {n};
end

if isempty(np)
    error('No point-matches found');
end
clear sID_all
disp('... concatenating point matches ...');
% concatenate
M = vertcat(M{:});
adj = vertcat(adj{:});
W   = vertcat(W{:});
np  = [np{:}]';

%%% additionally we return the vector of disconnected tiles that
%%% might need to be discarded
discard = setdiff(1:ntiles, unique(adj(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [M, W, adj, np] = filter_pm_local(M, W, adj, np, opts)
warning off;

if nargin < 2
    opts.NumRandomSamplingsMethod = 'Desired confidence';
    opts.MaximumRandomSamples = 3000;
    opts.DesiredConfidence = 99.5;
    opts.PixelDistanceThreshold = 1;
end

verbose = 0;

%% filter point matches using RANSAC
%warning('off', 'vision:obsolete:obsoleteFunctionality');
geoTransformEst = vision.GeometricTransformEstimator; % defaults to RANSAC
geoTransformEst.Method = 'Least Median of Squares';%'Random Sample Consensus (RANSAC)'; %
geoTransformEst.Transform = 'Nonreflective similarity';%'Affine'; % Valid values: 'Affine',
geoTransformEst.NumRandomSamplingsMethod = opts.NumRandomSamplingsMethod;% 'Desired confidence';
geoTransformEst.MaximumRandomSamples = opts.MaximumRandomSamples;%3000;
geoTransformEst.DesiredConfidence = opts.DesiredConfidence; %99.5;
geoTransformEst.PixelDistanceThreshold = opts.PixelDistanceThreshold; %0.01;


% M = pm.M;
% W = pm.W;
% adj = pm.adj;
% np  = pm.np;
del_ix = zeros(size(M,1),1);
parfor pmix = 1:size(M, 1)
    %warning('off', 'vision:obsolete:obsoleteFunctionality');
    m = M(pmix,:);
    if size(m{1},1)<2   % less than two points is useless
        del_ix(pmix) = pmix;
    else
        m1 = m{1};
        m2 = m{2};
        [tform_matrix, inlierIdx] = step(geoTransformEst, m{2}, m{1});
        m1 = m1(inlierIdx,:);
        m2 = m2(inlierIdx,:);
        M(pmix,:) = {m1,m2};
        w = W{pmix};
        if verbose > 1
            if sum(inlierIdx) < numel(inlierIdx)
                disp(['Eliminated: ' num2str(sum(~(inlierIdx))) ' points between ' mat2str(adj(pmix, :)) ' - remaining ' num2str(sum(inlierIdx))]);
            end
        end
        if sum(inlierIdx) == 0
            if verbose > 1
                disp(['Marked for deletion the link between ' mat2str(adj(pmix, :))])
            end
            del_ix(pmix) = pmix;
        end
        W{pmix} = w(inlierIdx);
        np(pmix) = length(W{pmix});
    end
end
del_ix(del_ix==0) = [];
M(del_ix,:) = [];
W(del_ix) = [];
adj(del_ix,:) = [];
np(del_ix) = [];


