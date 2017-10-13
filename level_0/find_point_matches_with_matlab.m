function [m12_2, m12_1, fh] = find_point_matches_with_matlab(im1,im2, method, visible, suppressError, opt)
if nargin<4, visible = true; end
if nargin<5, suppressError = false; end
if nargin<6, opt = []; end
if ~isfield(opt, 'verbose'), opt.verbose = 0; end
B = [];
T = [];
if ~isfield(opt, 'SURF_NumOctaves'), opt.SURF_NumOctaves = 1; end
if ~isfield(opt, 'SURF_NumScaleLevels'), opt.SURF_NumScaleLevels = 10; end
if ~isfield(opt, 'SURF_MetricThreshold'), opt.SURF_MetricThreshold = 100; end

if ~isfield(opt, 'MatchThreshold'), opt.MatchThreshold = 100; end
if ~isfield(opt, 'Unique'), opt.Unique = true; end
if ~isfield(opt, 'Method'), opt.Method = 'Exhaustive'; end
if ~isfield(opt, 'MaxRatio'), opt.MaxRatio = 0.65; end
if ~isfield(opt, 'Metric'), opt.Metric = 'SAD'; end

if opt.verbose
    disp('SURF and match parameters:');
    disp(opt);
end

[f1, vp1] = im_get_features(im2,method, opt);
[f2, vp2] = im_get_features(im1,method, opt);

[m12_1, m12_2] = im_pair_match_features_local(f1, vp1, f2, vp2, suppressError, opt);
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

function [m1, m2] = im_pair_match_features_local(f1, vp1, f2, vp2, suppressError, opt)
%
%%%% m1 and m2 are corresponding point group objects that were generated
%%%% without any assumption about the transformation, but were selected
%%%% based on the assumption of a specific transformation

index_pairs = matchFeatures(f1, f2,...
    'MatchThreshold', opt.MatchThreshold,...
    'Unique', opt.Unique, ...
    'Method', opt.Method, ...
    'MaxRatio', opt.MaxRatio,...
    'Metric', opt.Metric);

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