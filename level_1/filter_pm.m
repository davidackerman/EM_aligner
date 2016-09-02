function pm = filter_pm(pm, opts)
warning off;

if nargin<2
    opts.NumRandomSamplingsMethod = 'Desired confidence';
    opts.MaximumRandomSamples = 3000;
    opts.DesiredConfidence = 99.5;
    opts.PixelDistanceThreshold = 0.01;
end


warning off
% % %% filter point matches using RANSAC
% geoTransformEst = vision.GeometricTransformEstimator; % defaults to RANSAC
% geoTransformEst.Method = 'Random Sample Consensus (RANSAC)'; %'Least Median of Squares';%
% geoTransformEst.Transform = 'Nonreflective similarity';%'Affine';%
% geoTransformEst.NumRandomSamplingsMethod = 'Desired confidence';
% geoTransformEst.MaximumRandomSamples = 1000;
% geoTransformEst.DesiredConfidence = 99.99;
% geoTransformEst.PixelDistanceThreshold = 1.0;

% %% filter point matches using RANSAC
geoTransformEst = vision.GeometricTransformEstimator; % defaults to RANSAC
geoTransformEst.Method = 'Random Sample Consensus (RANSAC)'; %'Least Median of Squares';%
geoTransformEst.Transform = 'Nonreflective similarity';%'Affine';
geoTransformEst.NumRandomSamplingsMethod = opts.NumRandomSamplingsMethod;% 'Desired confidence';
geoTransformEst.MaximumRandomSamples = opts.MaximumRandomSamples;%3000;
geoTransformEst.DesiredConfidence = opts.DesiredConfidence; %99.5;
geoTransformEst.PixelDistanceThreshold = opts.PixelDistanceThreshold; %0.01;
warning on;

M= pm.M;
W = pm.W;
adj = pm.adj;
np  = pm.np;
del_ix = zeros(size(pm.M,1),1);
parfor pmix = 1:size(pm.M,1)
    m = M(pmix,:);
    m1 = m{1};
    m2 = m{2};
    [tform_matrix, inlierIdx] = step(geoTransformEst, m{2}, m{1});
    m1 = m1(inlierIdx,:);
    m2 = m2(inlierIdx,:);
    M(pmix,:) = {m1,m2};
    w = W{pmix};
%     if sum(inlierIdx)<numel(inlierIdx),
%         disp(['Eliminated: ' num2str(sum(~(inlierIdx))) ' points.']);
%     end
    if sum(inlierIdx)==0
        %del_ix = [del_ix;pmix];
        del_ix(pmix) = pmix;
    end
    W{pmix} = w(inlierIdx);
    np(pmix) = length(W{pmix});
end
warning on;
del_ix(del_ix==0) = [];
M(del_ix,:) = [];
W(del_ix) = [];
adj(del_ix,:) = [];
np(del_ix) = [];
pm.M = M;
pm.W = W;
pm.adj = adj;
pm.np = np;
