function [T, B] = register_image_pair(fixed,moving)

B = [];
T = [];
method = 'SURF';%'MSER_SURF';%'SURF';%'SURF'; %'BRISK' 'HARRIS' 'FAST''FAST_SURF';%
opt.SURF_NumOctaves = 4;
opt.SURF_NumScaleLevels = 16;
opt.SURF_MetricThreshold = 300;

[f1, vp1] = im_get_features(moving,method, opt);
[f2, vp2] = im_get_features(fixed,method, opt);

[m12_1, m12_2, tf12] = im_pair_match_features(f1, vp1, f2, vp2);
warning off;
thresh = 2;
maxIm = 255;
if ~isempty(m12_1) && ~isempty(m12_2)
    T = fitgeotrans([m12_1.Location],[m12_2.Location],'similarity');
    B = imwarp(moving, T,'FillValues', maxIm, 'OutputView',  imref2d(size(fixed)));
    %B = imwarp(im2, T,'FillValues', maxIm);
    
    %%%%%%%%% sosi
    figure; showMatchedFeatures(fixed, moving, m12_2, m12_1, 'montage');

end