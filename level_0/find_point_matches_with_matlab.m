function [m12_2, m12_1, fh] = find_point_matches_with_matlab(im1,im2, method, visible, suppressError)
if nargin<4, visible = true; end
if nargin<5, suppressError = false; end
B = [];
T = [];
%method = 'SURF';%'MSER_SURF';%'SURF';%'SURF'; %'BRISK' 'HARRIS' 'FAST''FAST_SURF';%
opt.SURF_NumOctaves = 1;
opt.SURF_NumScaleLevels = 10;
opt.SURF_MetricThreshold = 100;

[f1, vp1] = im_get_features(im2,method, opt);
[f2, vp2] = im_get_features(im1,method, opt);

[m12_1, m12_2] = im_pair_match_features_local(f1, vp1, f2, vp2, suppressError);
warning off;
thresh = 2;
maxIm = 255;

fh = figure('visible', visible);
if ~isempty(m12_1) % matching features found
    showMatchedFeatures(im1, im2, m12_2, m12_1, 'montage');
else
    imshow(imfuse(im1, im2, 'montage'));
end
end

function [m1, m2] = im_pair_match_features_local(f1, vp1, f2, vp2, suppressError)
%
%%%% m1 and m2 are corresponding point group objects that were generated
%%%% without any assumption about the transformation, but were selected
%%%% based on the assumption of a specific transformation

index_pairs = matchFeatures(f1, f2,...
    'MatchThreshold', 99,...
    'Unique', true, ...
    'Method', 'Exhaustive', ...
    'Metric', 'SAD');
if isempty(index_pairs)
    if suppressError
        warning on;
        warning('No matches found');
        warning off;
    else
        error('No matches found');
    end
end

m1  = vp1(index_pairs(:,1));
m2  = vp2(index_pairs(:,2));
end