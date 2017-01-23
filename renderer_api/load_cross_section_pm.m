function [M, adj, W, np, xnp] = load_cross_section_pm(pm, sID1, sID2, map_id, min_points, max_points, wopts, fac)
% returns point-matches (from pm) between two sections with sID1 and sID2.
% sID{ix-(n-1)} is sID2, where n>=2 is the number of sections we are removed from sID1 (based on its z
% value). sID2<sID1
M = [];
adj = [];
W = [];
np = [];
xnp = 0;
count = 1;
for six1 = 1:numel(sID1)                % loop over reacquires if present. Usually numel(sID1) is equal to 1
    for six2 = 1:numel(sID2)            % loop over reacquires if present. Usually numel(sID2) is equal to 1
        
        j = get_pms_cross_layer(pm, sID2{six2}, sID1{six1}, wopts);

        for jix = 1:numel(j)
            if size(j(jix).matches.p',1)>=min_points
                % check that point-match is between tiles/canvases that exist in our working set
                % i.e. only ids that exist in map_id will be used
                if isKey(map_id, j(jix).pId) && isKey(map_id, j(jix).qId)
                    % obtain map ids
                    midp = map_id(j(jix).pId);
                    midq = map_id(j(jix).qId);
                    if numel(j(jix).matches.p(1,:))>max_points
                        indx = randi(numel(j(jix).matches.p(1,:))-1, max_points,1);
                        if midp<midq
                            M{count,1}   = j(jix).matches.p(1:2,indx)';
                            M{count,2}   = j(jix).matches.q(1:2,indx)';
                            adj(count,:) = [midp midq];
                        else
                            M{count,1}   = j(jix).matches.q(1:2,indx)';
                            M{count,2}   = j(jix).matches.p(1:2, indx)';
                            adj(count,:) = [midq midp];
                        end
                        W{count,1}   = fac * j(jix).matches.w(indx)';         % relative weights of point matches within this group
                        np(count)  = max_points;    %  we are recording the number of point matches between those two tiles
                    else
                        if midp<midq
                            M{count,1}   = j(jix).matches.p';
                            M{count,2}   = j(jix).matches.q';
                            adj(count,:) = [midp midq];
                        else
                            M{count,1}   = j(jix).matches.q';
                            M{count,2}   = j(jix).matches.p';
                            adj(count,:) = [midq midp];
                        end
                        W{count,1}   = fac * j(jix).matches.w';         % relative weights of point matches within this group
                        np(count)  = size(j(jix).matches.p',1);    %  we are recording the number of point matches between those two tiles
                    end
                    count = count + 1;
                    xnp = xnp +  1;
                end
            end
        end
    end
end