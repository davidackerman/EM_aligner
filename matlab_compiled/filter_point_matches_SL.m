function filter_point_matches_SL(fn)
sl = loadjson(fileread(fn));


%% filter point matches using RANSAC
%warning('off', 'vision:obsolete:obsoleteFunctionality');
geoTransformEst = vision.GeometricTransformEstimator; % defaults to RANSAC
geoTransformEst.Method = 'Random Sample Consensus (RANSAC)';%'Least Median of Squares';%'Random Sample Consensus (RANSAC)'; %
geoTransformEst.Transform = 'Nonreflective similarity';%'Affine'; % Valid values: 'Affine',
geoTransformEst.NumRandomSamplingsMethod = sl.pmopts.NumRandomSamplingsMethod;% 'Desired confidence';
geoTransformEst.MaximumRandomSamples = sl.pmopts.MaximumRandomSamples;%3000;
geoTransformEst.DesiredConfidence = sl.pmopts.DesiredConfidence; %99.5;
geoTransformEst.PixelDistanceThreshold = sl.pmopts.PixelDistanceThreshold; %0.01;


% M = pm.M;
% W = pm.W;
% adj = pm.adj;
% np  = pm.np;
canvasPairs = sl.canvasPairs';
del_ix = zeros(size(canvasPairs,1),1);
for pmix = 1:size(canvasPairs, 1)
    matches_p = canvasPairs{pmix}.matches.p';
    matches_q = canvasPairs{pmix}.matches.q';
    if size(matches_p,1)<2   % less than two points is useless
        del_ix(pmix) = pmix;
    else
        [tform_matrix, inlierIdx] = step(geoTransformEst, matches_q, matches_p);
        canvasPairs{pmix,:}.matches.p = matches_p(inlierIdx,:)';
        canvasPairs{pmix,:}.matches.q = matches_q(inlierIdx,:)';
%        if verbose > 1
%            if sum(inlierIdx) < numel(inlierIdx)
%                disp(['Eliminated: ' num2str(sum(~(inlierIdx))) ' points between ' mat2str(adj(pmix, :)) ' - remaining ' num2str(sum(inlierIdx))]);
%            end
%        end
        if sum(inlierIdx) == 0
%            if verbose > 1
%                disp(['Marked for deletion the link between ' mat2str(adj(pmix, :))])
%            end
            del_ix(pmix) = pmix;
        end
        canvasPairs{pmix,:}.matches.w = canvasPairs{pmix,:}.matches.w(inlierIdx);
%       np(pmix) = length(W{pmix});
    end
end
del_ix(del_ix==0) = [];
canvasPairs(del_ix) = [];
json_string = savejson('',canvasPairs);
fprintf('canvasPairs: \n');
fprintf(json_string);
%W(del_ix) = [];
%adj(del_ix,:) = [];
%np(del_ix) = [];