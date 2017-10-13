%% tile information

% render collection

rc.baseURL = 'http://tem-services:8080/render-ws/v1';

rc.owner = 'hessh';

rc.project = 'fly_em_pre_iso';

rc.stack = 'column_30_acquire';

%tile ids

tile_1_id = '17-02-11_014625_0-0-0.23279.0';

tile_2_id = '17-02-11_014625_0-0-1.23279.0';

%% url parameters used for rendering tiles, and the defaults that will be used if not provided

% DEFAULT VALUES:

%   normalizeForMatching    : true

%   renderWithFilter        : true

%   renderWithoutMask       : true

%   fullScaleWidth          : unset (null)

%   fullScaleHeight         : unset (null)

 

% url options struct

url_options.renderWithFilter = 'false';
url_options.renderWithoutMask = 'true';

%% SIFT parameters and the defaults that will be used if not provided

%     SIFTfdSize

%        SIFT feature descriptor size: how many samples per row and column

%        Default: 89

%     SIFTmaxScale

%        SIFT maximum scale: minSize * minScale < size < maxSize * maxScale

%        Default: 0.85

%     SIFTminScale

%        SIFT minimum scale: minSize * minScale < size < maxSize * maxScale

%        Default: 0.5

%     SIFTsteps

%        SIFT steps per scale octave

%        Default: 3

%     fillWithNoise

%        Fill each canvas image with noise before rendering to improve point match

%        derivation

%        Default: true

%     matchIterations

%        Match filter iterations

%        Default: 1000

%     matchMaxEpsilon

%        Minimal allowed transfer error for match filtering

%        Default: 20.0

%     matchMaxNumInliers

%        Maximum number of inliers for match filtering

%     matchMaxTrust

%        Reject match candidates with a cost larger than maxTrust * median cost

%        Default: 3.0

%     matchMinInlierRatio

%        Minimal ratio of inliers to candidates for match filtering

%        Default: 0.0

%     matchMinNumInliers

%        Minimal absolute number of inliers for match filtering

%        Default: 4

%     matchModelType

%        Type of model for match filtering

%        Default: AFFINE

%        Possible Values: [TRANSLATION, RIGID, SIMILARITY, AFFINE]

%     matchRod

%        Ratio of distances for matches

%        Default: 0.92

%     renderScale

%        Render canvases at this scale

%        Default: 1.0

%     outputDirectory

%        Parent directory in which subdirectories will be created to store

%        images and point-match results from SIFT.

%        NO DEFAULT VALUE

% SIFT options struct. To test multiple values for a single parameter,

% set it using an array

SIFT_options = [];

% SIFT_options.SIFTfdSize = 8;
% 
% SIFT_options.SIFTminScale = 0.58;
% 
% SIFT_options.renderScale = 0.6;
% 
% SIFT_options.SIFTmaxScale = 0.6;
% 
% SIFT_options.matchModelType = 'TRANSLATION';
% 
% SIFT_options.matchMaxEpsilon = 5.0;
% 
% SIFT_options.matchMinNumInliers = 4;
% 
% SIFT_options.matchMaxNumInliers = 50;
% 
% SIFT_options.outputDirectory = [pwd '/SIFT_output/'];

 

%% SURF parameters, currently only renderScale is used
SURF_options.verbose = 1;
SURF_options.renderScale = 1.0;
SURF_options.numPixels = 550;   % cut off window size for tiles to be compared

%%%% SURF-specific parameters
SURF_options.SURF_NumOctaves = 2;
SURF_options.SURF_NumScaleLevels = 6;
SURF_MetricThreshold = 1;    % decrease to return more features

%%%% match parameters
SURF_options.MatchThreshold = 100;   % increase to return more matches (max value is 100)
SURF_options.Unique = true;
SURF_options.Method = 'Exhaustive';
SURF_options.MaxRatio = 0.7;
SURF_options.Metric = 'SAD';

%% Set the output directory for storing the final images, as well as the desired visibility of the figures, and run the code
clc;
result_output_directory = [pwd '/point_match_parameterization_result_output'];

figure_visibility = true;

% if SIFT_options or SURF_options is empty, then it is skipped.

point_match_optimization(rc, tile_1_id, tile_2_id, SIFT_options, SURF_options, url_options, result_output_directory, figure_visibility);



%%


















































