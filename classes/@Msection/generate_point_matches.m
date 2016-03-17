function [obj, M, adj] = generate_point_matches(obj, thresh, filter)
% calculates all features for all tiles and generates point matches
% structure, with point locations, adjacency and weights
% SOSI: parallel point-matches generation is needed
% SOSI: does not handle orphan tiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2, thresh = 12;end;
if nargin<3, filter = 'true';end

if isempty(obj.A), 
    warning('No adjacency matrix found, calculating');
    obj = update_XY(obj);
    obj = update_adjacency(obj);
end
if size(obj.A,1)~=numel(obj.tiles)
    warning('Invalid or outdated adjacency matrix, recalculating');
    obj = update_XY(obj);
    obj  = update_adjacency(obj);
end

[obj] = calculate_tile_features(obj, filter);

[r, c] = ind2sub(size(obj.A), find(obj.A));  % neighbors are determined by the adjacency matrix

%% % distributed tile-pair matching
mL2_tiles = obj.tiles;
M = cell(numel(r),2);
adj = zeros(numel(r),2);
W = cell(numel(r),1);
np = zeros(numel(r),1);
delpix = zeros(numel(r),1, 'uint32');
%count = 1;
disp('Calculating point matches .... ');
%parfor_progress(numel(r));
parfor pix = 1: numel(r)
    %disp(['Point matching: ' num2str(pix) ' of ' num2str(numel(r))]);
    f1 = mL2_tiles(r(pix)).features;
    f2 = mL2_tiles(c(pix)).features;
    vp1 = mL2_tiles(r(pix)).validPoints;
    vp2 = mL2_tiles(c(pix)).validPoints;
    [m1, m2, ~]  = im_pair_match_features(f1, vp1, f2, vp2);
    
    M(pix,:) = {[m1.Location],[m2.Location]};
    adj(pix,:) = [r(pix) c(pix)];
    W(pix) = {[ones(size(m1.Location,1),1) * 1/(1+ abs(mL2_tiles(r(pix)).z-mL2_tiles(c(pix)).z))]};
    np(pix)  = size(m1.Location,1);
    %%%% mark for removal point-matches that don't have enough point pairs
    if size(m1.Location,1)<thresh
        delpix(pix) = 1;
    end
%    parfor_progress;
end
%parfor_progress(0);
if isempty(M), error('No matches found');end;
disp('Done!');
delpix = logical(delpix);
M(delpix,:) = [];
adj(delpix,:) = [];
W(delpix) = [];
np(delpix) = [];

obj.pm.M = M;
obj.pm.adj = adj;
obj.pm.W = W;
obj.pm.np = np;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %% generates lists of matching points M for each pair of adjacent tiles and returns the pairs 2-vector adj
% % %%%% needs to be distributed
% % verbose = 1;
% % %%% Note: Deprecated ----> use generate_features_point_matches instead
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %% %% get SURF features for each image
% % method = 'SURF';
% % f = {};
% % vp = {};
% % disp('Getting features....');
% % for tix = 1:numel(obj.tiles)
% %     disp(tix);
% %     [f{tix}, vp{tix}] = get_features(obj.tiles(tix).get_image(),method);% % extract matching point pairs using Matlab's SURF
% % end
% % disp('Done!');
% % 
% % %% obtain matching point pairs for each set of overlaping tiles
% % 
% % [r, c] = ind2sub(size(obj.A), find(obj.A(:)));
% % pcount = 0;
% % M = {};
% % adj = [];
% % disp('Generating point matches...');
% % for ix = 1:numel(r)
% %     if verbose, disp([num2str(ix) ' of ' num2str(numel(r))]);end
% % 
% %     [m1, m2]= im_pair_match_features(f{r(ix)}, vp{r(ix)}, f{c(ix)}, vp{c(ix)});
% %     
% %         pcount = pcount + 1;
% %         M{pcount, 1} = m1;
% %         M{pcount, 2} = m2;
% %         adj = [adj; r(ix) c(ix)];
% %    
% % end
% % disp('Done!');































