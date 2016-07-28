function pm = filter_pm(pm)
warning off;
% %% filter point matches using RANSAC
geoTransformEst = vision.GeometricTransformEstimator; % defaults to RANSAC
geoTransformEst.Method = 'Random Sample Consensus (RANSAC)'; %'Least Median of Squares';%
geoTransformEst.Transform = 'Nonreflective similarity';%'Affine';%
geoTransformEst.NumRandomSamplingsMethod = 'Desired confidence';
geoTransformEst.MaximumRandomSamples = 1000;
geoTransformEst.DesiredConfidence = 99.99;

geoTransformEst.PixelDistanceThreshold = 1.0;

M = pm.M;
W = pm.W;
adj = pm.adj;
np  = pm.np;
del_ix = [];
for pmix = 1:size(pm.M,1)

    
    m1 = M{pmix,1};
    m2 = M{pmix,2};
    [tform_matrix, inlierIdx] = step(geoTransformEst, m2, m1);
    m1 = m1(inlierIdx,:);
    m2 = m2(inlierIdx,:);
    M{pmix,1} = m1;
    M{pmix,2} = m2;
    w = W{pmix};
%     if sum(inlierIdx)<numel(inlierIdx),
%         disp(['Eliminated: ' num2str(sum(~(inlierIdx))) ' points.']);
%     end
    if sum(inlierIdx)==0
        del_ix = [del_ix;pmix];
    end
    W{pmix} = w(inlierIdx);
    np(pmix) = length(W{pmix});
end
warning on;
M(del_ix,:) = [];
W(del_ix) = [];
adj(del_ix,:) = [];
np(del_ix) = [];
pm.M = M;
pm.W = W;
pm.adj = adj;
pm.np = np;
