function [L, tIds, PM, pm_mx, sectionId, z] = load_point_matches(nfirst, ...
    nlast, rc, pm, nbr, min_points, xs_weight, max_points)
% Input: nfirst and nlast are zvalue of sections in rc
%        rc and pm are structs with specifications for accessing
%        collections of tile-specs and point-matches (that are related)
%        see exmple below.
%        (optional)
%        nbr: number of neighboring sections to consider
%        min_points: minimum number of points between two tiles
%        xs_weight: weight factor for cross-section point matches
% Output: Msection object L with field pm (which is a struct with fields M, adj, W and np), and tileIds
%         M: is a cell array of size npx2, e.g. a set of point matches is given by M{1,1} for xy of
%         the first set of points and M{1,2} for xy of the matching points
%
% requires a starting value for section zvalues (nfirst) and the zvalue of
% the last section (nlast). nfirst and nlast can be the same value
% First a sectionID list is created in the order of zvalues (most of those
% are already ordered, but we need to be sure). This is why we need rc.
% Second all point matches within the sections and across sections will be
% downloaded from the pm database.
%
% % Example input:
% % rc.stack = 'v12_align';
% % rc.owner='flyTEM';
% % rc.project='FAFB00';
% % rc.server='http://tem-services.int.janelia.org:8080/render-ws/v1';
% %
% % pm.server = 'http://tem-services.int.janelia.org:8080/render-ws/v1';
% % pm.owner  = 'flyTEM';
% % pm.match_collection = 'v12_dmesh';
% % nbr = 3;
% % min_points = 5;
% % xs_weight = 1;
%
% Author: Khaled Khairy. Janelia Research Campus 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<5, nbr = 4;end  % number of neighbors to check
if nargin<6, min_points = 0;end
if nargin<7, xs_weight = 1;end
if nargin<8, max_points = inf;end
% %% get the list of zvalues and section ids within the z range between nfirst and nlast (inclusive)
% urlChar = sprintf('%s/owner/%s/project/%s/stack/%s/sectionData', ...
%     rc.baseURL, rc.owner, rc.project, rc.stack);
% js = webread(urlChar);
% sectionId = {js(:).sectionId};
% [z, ia]   = sort(([js(:).z]));
% sectionId = sectionId(ia);
%
%
% indx = find(z>=nfirst & z<=nlast);
% sectionId = sectionId(indx);% determine the sectionId list we will work with
% z         = z(indx);        % determine the zvalues (this is also the spatial order)
%
% % we need unique values, and we need to know how many sectionId's correspond to each unique z value
% % usually it is one, but sometimes we have hi/lo dose or other regions
% [zu, ia, ic] = unique(z);
% count = 1;
% sID = {};
% for zix = 1:numel(zu)
%     ns(zix) =  numel(find(ic==zix));
%     vec = {};
%     for six = 1:ns(zix)
%         vec{six} = sectionId{count};
%         sID{zix} = vec;
%         count = count + 1;
%     end
% end
[zu, sID, sectionId, z, ns] = get_section_ids(rc, nfirst, nlast);
%% get a list of all tiles for those sections
options = weboptions;
options.Timeout = 20;
clear t;
tilecount = [];
parfor ix = 1:numel(zu)
    urlChar = sprintf('%s/owner/%s/project/%s/stack/%s/z/%d/tile-specs', ...
        rc.baseURL, rc.owner, rc.project, rc.stack, zu(ix));
    j = webread(urlChar, options);
    jt = tile;
    tilecount(ix) = numel(j);
    for jix = 1:numel(j)
        jt(jix) = tile(j(jix));
        jt(jix).z = zu(ix);
        
    end
    t(ix).jt = jt;
end
ntiles = sum(tilecount);
% concatenate all tile ids
tIds = cell(ntiles,1);
tiles(ntiles) = tile;
cnt = 1;
for ix = 1:numel(zu)
    tIds(cnt:cnt+tilecount(ix)-1) = {t(ix).jt.renderer_id};
    tiles(cnt:cnt+tilecount(ix)-1) = t(ix).jt;
    cnt = cnt + tilecount(ix);
    %     tIds = [tIds {t(ix).jt.renderer_id}];
    %     tiles = [tiles t(ix).jt];
end

% loop over tiles to set tile id
for ix = 1:numel(tiles)
    tiles(ix).id = ix;
    tiles(ix).owner = rc.owner;
    tiles(ix).project = rc.project;
    tiles(ix).stack = rc.stack;
end
L = Msection(tiles);
%L = update_tile_sources(L, rc);

% %%%%%%%%%%%% check consistency
% % check that all renderer_ids in L are unique (this means that all tiles are unique)
% % if not then flag which ones are not
% rids = {L.tiles(:).renderer_id};
% [un idx_last idx] = unique(rids);
% uqindx = accumarray(idx(:),(1:length(idx))',[],@(x) {sort(x)});
% for ix = 1:numel(uqindx)
%     if numel(uqindx{ix})>1,
%         id = L.tiles(uqindx{ix}(1)).renderer_id;
%         disp([num2str(numel(uqindx{ix})) ' copies of id: ' id ' found.']);
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% map renderer ids to their index position in L
clear count_vec ;
clear id_vec ;
parfor ix = 1:numel(tIds)
    count_vec(ix)= {ix};
    id_vec(ix) = tIds(ix);%tIds{ix};
end
map_id = containers.Map(id_vec, count_vec);

%% get point matches for each section (montage)
PM.M = {};
PM.adj = [];
PM.W = {};
PM.np = [];
n1 = [];
parfor ix = 1:numel(ns)
    %disp(ix);
    count = 1;
    n1(ix) = 0;
    for six = 1:ns(ix)
        %disp([six count]);
        urlChar = sprintf('%s/owner/%s/matchCollection/%s/group/%s/matchesWithinGroup', ...
            pm.server, pm.owner, pm.match_collection, sID{ix}{six});
        try
            jj = webread(urlChar, options);
        catch err_fetch_pm
            kk_disp_err(err_fetch_pm)
            pause(1);
            jj = webread(urlChar,options); % try again
        end
        n1(ix) = n1(ix) + numel(jj);
        for jix = 1:numel(jj)
            if size(jj(jix).matches.p',1)>=min_points
                if isKey(map_id, jj(jix).pId) && isKey(map_id, jj(jix).qId)
                    if numel(jj(jix).matches.p(1,:))>max_points
                        indx = randi(numel(jj(jix).matches.p(1,:))-1, max_points,1);
                        PM(ix).M{count,1}   = [jj(jix).matches.p(1:2,indx)]';
                        PM(ix).M{count,2}   = [jj(jix).matches.q(1:2,indx)]';
                        PM(ix).adj(count,:) = [map_id(jj(jix).pId) map_id(jj(jix).qId)];
                        PM(ix).W{count,1}     = jj(jix).matches.w(indx)';         % relative weights of point matches within this group
                        PM(ix).np(count)    = max_points;
                    else
                        PM(ix).M{count,1}   = [jj(jix).matches.p]';
                        PM(ix).M{count,2}   = [jj(jix).matches.q]';
                        PM(ix).adj(count,:) = [map_id(jj(jix).pId) map_id(jj(jix).qId)];
                        PM(ix).W{count,1}     = jj(jix).matches.w';         % relative weights of point matches within this group
                        PM(ix).np(count)    = size(jj(jix).matches.p',1);
                    end
                    
                    count = count + 1;
                end
            end
        end
    end
end


%% obtain cross-section point-matches
%disp('Obtaining cross-layer point_matches');
xPM = {};
n   = {};
parfor pmix = 1:nbr
    [xPM{pmix}, n{pmix}] = get_cross_section_pm(pmix+1, ...
        pm, sID, map_id, min_points, xs_weight, max_points);%% get point matches to immediate neighbor
end
%disp('Generating final point match struct');
%% generate final M, adj, W and np
M   = [];
adj = [];
W   = [];
np = [];
for ix = 1:numel(zu)   % loop over sections
    %disp(ix);
    M = [M;PM(ix).M];
    adj = [adj;PM(ix).adj];
    W   = [W;PM(ix).W];
    np   = [np;PM(ix).np(:)];
    for nix = 1:nbr   % loop over neighboring sections
        %disp([ix nix]);
        if  ~(numel(xPM{nix})==1 && isempty(xPM{nix}.M))
            if numel(xPM{nix})>=ix
                if ~isempty(xPM{nix}(ix).M)
                    M = [M;xPM{nix}(ix).M];
                    adj = [adj;xPM{nix}(ix).adj];
                    W = [W;xPM{nix}(ix).W];
                    np = [np;xPM{nix}(ix).np(:)];
                end
            end
        end
    end
end

%% place point match information into Msection object L
L.pm.M = M;
L.pm.adj = adj;
L.pm.W = W;
L.pm.np = np;
% %% tiles that are not connected to any other tiles must be eliminated
% adjids = {L.tiles(adj(:)).renderer_id};
%% %%%%%%%%%%%%%%%%%%%%%%%%% diagnostics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% section adjacency point-match count matrix
pm_mx = diag(n1);
for nbrix = 1:nbr
    if ~isempty(n{nbrix})
        if ~(sum(n{nbrix}(:,1))==0),
            pm_mx = pm_mx + diag(n{nbrix}(:,1),nbrix);
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%% Check consistency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(sum(pm_mx)==0), disp('Warning: defective pm connectivity matrix');end
% L.pm.adj are unique: report error if not
% --------- sosi: invesitgate how duplicates could arise in the first place
[bb, indx] = unique(L.pm.adj,'rows');
if ~(size(bb,1)==size(L.pm.adj,1))
    error('Rows in L.pm.adj should be unique');
end
%%%%%%%%%%%%%%%%%%%%%%%%
function [xPM, np_vec] = get_cross_section_pm(n, pm, sID, map_id, min_points, xs_weight, max_points)
% assumes sectionId contains sections sorted by z
% small xs_weight means less weight for cross-layer point-match
if nargin<6, xs_weight = 1;end
np_vec = [];
xPM.M = {};
xPM.adj = [];
xPM.W = {};
xPM.np = [];
options = weboptions;
options.Timeout = 20;
sec_ix = 1;
if numel(sID)>=n
    np_vec = zeros(numel(n:numel(sID)),1);
    fac = (1/n * xs_weight);   % small fac ==> less weight for cross-layer point matches
    
    for ix = n:numel(sID)
        np_vec(sec_ix) = 0;
        count = 1;
        for six1 = 1:numel(sID{ix})
            for six2 = 1:numel(sID{ix-(n-1)})
                urlChar = sprintf('%s/owner/%s/matchCollection/%s/group/%s/matchesWith/%s', ...
                    pm.server, pm.owner, pm.match_collection, sID{ix-(n-1)}{six2}, sID{ix}{six1});
                j = webread(urlChar, options);
                
                %         np_vec(ix-(n-1)) = numel(j);
                for jix = 1:numel(j)
                    if size(j(jix).matches.p',1)>=min_points
                        if isKey(map_id, j(jix).pId) && isKey(map_id, j(jix).qId)
                            
                            % %                             % sosi
                            % %                             if str2double(sID{ix}{six1})==10
                            % %                                 disp([map_id(j(jix).pId)<map_id(j(jix).qId) ...
                            % %                                      n (ix-(n-1)) ix map_id(j(jix).pId) map_id(j(jix).qId)]);
                            % %                             end
                            
                            % obtain map ids
                            midp = map_id(j(jix).pId);
                            midq = map_id(j(jix).qId);
                            
                            if numel(j(jix).matches.p(1,:))>max_points
                                indx = randi(numel(j(jix).matches.p(1,:))-1, max_points,1);
                                if midp<midq
                                    xPM(ix-(n-1)).M{count,1}   = j(jix).matches.p(1:2,indx)';
                                    xPM(ix-(n-1)).M{count,2}   = j(jix).matches.q(1:2,indx)';
                                    xPM(ix-(n-1)).adj(count,:) = [midp midq];
                                else
                                    xPM(ix-(n-1)).M{count,1}   = j(jix).matches.q(1:2,indx)';
                                    xPM(ix-(n-1)).M{count,2}   = j(jix).matches.p(1:2, indx)';
                                    xPM(ix-(n-1)).adj(count,:) = [midq midp];
                                end
                                xPM(ix-(n-1)).W{count,1}   = fac * j(jix).matches.w(indx)';         % relative weights of point matches within this group
                                xPM(ix-(n-1)).np(count)  = max_points;    %  we are recording the number of point matches between those two tiles
                                
                            else
                                if midp<midq
                                    xPM(ix-(n-1)).M{count,1}   = j(jix).matches.p';
                                    xPM(ix-(n-1)).M{count,2}   = j(jix).matches.q';
                                    xPM(ix-(n-1)).adj(count,:) = [midp midq];
                                else
                                    xPM(ix-(n-1)).M{count,1}   = j(jix).matches.q';
                                    xPM(ix-(n-1)).M{count,2}   = j(jix).matches.p';
                                    xPM(ix-(n-1)).adj(count,:) = [midq midp];
                                end
                                xPM(ix-(n-1)).W{count,1}   = fac * j(jix).matches.w';         % relative weights of point matches within this group
                                xPM(ix-(n-1)).np(count)  = size(j(jix).matches.p',1);    %  we are recording the number of point matches between those two tiles
                                
                                
                            end
                            
                            count = count + 1;
                            np_vec(sec_ix) = np_vec(sec_ix) +  1;
                            
                        end
                    end
                end
                %pause(0.2);
            end
        end
        sec_ix = sec_ix + 1;
    end
end










































