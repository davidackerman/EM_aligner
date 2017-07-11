function [x, y, err] = local_to_world(rc, tileId, x, y)
% coordinate mapping of a point from local coordiantes xyz to rc

err = 0;
% map from tile local to rc world and return new coordinate
url = sprintf('%s/owner/%s/project/%s/stack/%s/tile/%s/local-to-world-coordinates/%s,%s',...
    rc.baseURL, rc.owner, rc.project, rc.stack, tileId, x, y);

U = matlab.net.URI(url);
try
jj2 = webread(char(U));

x = jj2.world(1);
y = jj2.world(2);
catch err_read_point
    x = [];
    y = [];
    err = 1;
end