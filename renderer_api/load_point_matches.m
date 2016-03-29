function [L, tIds, PM] = load_point_matches(nfirst, nlast, rc, pm)
% Input: nfirst and nlast are zvalue of sections in rc
%        rc and pm are structs with specifications for accessing
%        collections of tile-specs and point-matches (that are related)
% returns Msection object L with field pm (which is a struct with fields M, adj, W and np), and tileIds
% requires a starting value for section zvalues (nfirst) and the zvalue of
% the last section required (nlast). nfirst and nlast can be the same value
% First a sectionID list is created in the order of zvalues (most of those
% are already ordered, but we need to be sure). This is why we need rc.
% Second all point matches within the sections and across sections will be
% downloaded from the pm database.
% % Example rc and pm
% % rc.stack = 'v9_acquire_LC_merged_2';
% % rc.owner='flyTEM';
% % rc.project='FAFB00';
% % rc.server='http://tem-services.int.janelia.org:8080/render-ws/v1';
% % 
% % pm.server = 'http://tem-services.int.janelia.org:8080/render-ws/v1';
% % pm.owner  = 'flyTEM';
% % pm.match_collection = 'v9_1';

%% get the list of zvalues and section ids within the z range between nfirst and nlast (inclusive)
urlChar = sprintf('%s/owner/%s/project/%s/stack/%s/sectionData', ...
     rc.baseURL, rc.owner, rc.project, rc.stack);
 j = webread(urlChar);
 sectionId = {j(:).sectionId};
 z         = [j(:).z];
 indx = find(z>=nfirst & z<=nlast);
 
 sectionId = sectionId(indx);% determine the sectionId list we will work with
 z         = z(indx);        % determine the zvalues (this is also the spatial order)
 
 %% get a list of all tiles for those sections
 
 
%  % <sosi---- > this is how it should be done in the future --
%  instantiation using rc and z has already been implemented, see
%  Msection.m
%  parfor ix = 1:numel(z)
%      L(ix) = Msection(rc,z(ix));
%  end
%  L = concatenate_tiles(L); % concatenate all sections
%  %%%% sosi />


options = weboptions;
options.Timeout = 20;
clear t;
parfor ix = 1:numel(z)
     urlChar = sprintf('%s/owner/%s/project/%s/stack/%s/z/%d/tile-specs', ...
              rc.baseURL, rc.owner, rc.project, rc.stack, z(ix));
      j = webread(urlChar, options);
      jt = tile;
      for jix = 1:numel(j)
          jt(jix) = tile(j(jix));
          jt(jix).z = z(ix);
      end
      t(ix).jt = jt;
 end
 
 % concatenate all tile ids
 tIds = {};
 tiles = [];
 for ix = 1:numel(z)
     tIds = [tIds {t(ix).jt.renderer_id}];
     tiles = [tiles t(ix).jt];
 end
 
  % loop over tiles to set tile id
 parfor ix = 1:numel(tiles)
     tiles(ix).id = ix;
 end
 L = Msection(tiles);
 L = update_tile_sources(L, rc);
 
 
 % map renderer ids to their index position in L
clear count_vec ;
clear id_vec ;
parfor ix = 1:numel(tIds)
    count_vec(ix)= {ix};
    id_vec(ix) = tIds(ix);%tIds{ix};
end
map_id = containers.Map(id_vec, count_vec);

%% get point matches for each section
PM.M = {};
PM.adj = [];
PM.W = {};
PM.np = [];
min_points = 5;    
parfor ix = 1:numel(z)
    urlChar = sprintf('%s/owner/%s/matchCollection/%s/group/%s/matchesWithinGroup', ...
               pm.server, pm.owner, pm.match_collection, sectionId{ix});
           try
               j = webread(urlChar, options);
           catch err_ip_address
               pause(1);
               j = webread(urlChar,options);
           end
    count = 1;
    for jix = 1:numel(j)
        if size(j(jix).matches.p',1)>=min_points
            if isKey(map_id, j(jix).pId) && isKey(map_id, j(jix).qId)
                PM(ix).M{count,1}   = j(jix).matches.p';
                PM(ix).M{count,2}   = j(jix).matches.q';
                PM(ix).adj(count,:) = [map_id(j(jix).pId) map_id(j(jix).qId)];
                PM(ix).W{count,1}     = j(jix).matches.w';         % relative weights of point matches within this group
                PM(ix).np(count)    = size(j(jix).matches.p',1);
                count = count + 1;
            end
        end
    end
end

[xPM2] = get_cross_section_pm(2, pm, sectionId, map_id, min_points);%% get point matches to immediate neighbor
[xPM3] = get_cross_section_pm(3, pm, sectionId, map_id, min_points);%% get point matches one section removed from neighbor
%[xPM4] = get_cross_section_pm(4, pm, sectionId, map_id, min_points);%% get point matches two sectionss removed from neighbor


%% concatenate PM and xPM in the order that will be used for filling A in the solution
clear M adj W np;
for ix = 1:numel(z)
    if ix ==1
        M   = PM(ix).M;
        adj = PM(ix).adj;
        W   = PM(ix).W;
        np = PM(ix).np(:);
    elseif ix ==2
        M   = [M;   xPM2(ix-1).M;PM(ix).M];
        adj = [adj; xPM2(ix-1).adj;PM(ix).adj];
        W   = [W;   xPM2(ix-1).W;PM(ix).W];
        np =  [np;  xPM2(ix-1).np(:);PM(ix).np(:)];
    elseif numel(xPM3)>=(ix-2)
        M   = [M;   xPM3(ix-2).M;    xPM2(ix-1).M;PM(ix).M];
        adj = [adj; xPM3(ix-2).adj;  xPM2(ix-1).adj;PM(ix).adj];
        W   = [W;   xPM3(ix-2).W;    xPM2(ix-1).W;PM(ix).W];
        np =  [np;  xPM3(ix-2).np(:);xPM2(ix-1).np(:);PM(ix).np(:)];
    else
        M   = [M;   xPM2(ix-1).M;PM(ix).M];
        adj = [adj; xPM2(ix-1).adj;PM(ix).adj];
        W   = [W;   xPM2(ix-1).W;PM(ix).W];
        np =  [np;  xPM2(ix-1).np(:);PM(ix).np(:)];
    end
end


L.pm.M = M;
L.pm.adj = adj;
L.pm.W = W;
L.pm.np = np;


%%%%%%%%%%%%%%%%%%%%%%%%
function [xPM] = get_cross_section_pm(n, pm, sectionId, map_id, min_points)

xPM.M = {};
xPM.adj = [];
xPM.W = {};
xPM.np = [];     
options = weboptions;
options.Timeout = 20;
if numel(sectionId)>=n
    fac = (1/n)^2;
    for ix = n:numel(sectionId)
        urlChar = sprintf('%s/owner/%s/matchCollection/%s/group/%s/matchesWith/%s', ...
            pm.server, pm.owner, pm.match_collection, sectionId{ix-(n-1)}, sectionId{ix});
        j = webread(urlChar, options);
        count = 1;
        for jix = 1:numel(j)
            if size(j(jix).matches.p',1)>=min_points
                if isKey(map_id, j(jix).pId) && isKey(map_id, j(jix).qId)
                    xPM(ix-(n-1)).M{count,1}   = j(jix).matches.p';
                    xPM(ix-(n-1)).M{count,2}   = j(jix).matches.q';
                    xPM(ix-(n-1)).adj(count,:) = [map_id(j(jix).pId) map_id(j(jix).qId)];
                    xPM(ix-(n-1)).W{count,1}   = fac * j(jix).matches.w';         % relative weights of point matches within this group
                    xPM(ix-(n-1)).np(count)  = size(j(jix).matches.p',1);    %  we are recording the number of point matches between those two tiles
                    count = count + 1;
                end
            end
        end
        pause(1);
    end
end










































