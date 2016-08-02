function [x, y, tids, L1] = get_tile_centers(rc, z)
%%%% fast way to get tile ids, centers and minimal Msection object
tids = {};
urlChar = sprintf('%s/owner/%s/project/%s/stack/%s/z/%.1f/tile-specs', ...
    rc.baseURL, rc.owner, rc.project, rc.stack,z );
j = webread(urlChar);
jt1 = tile;
sectionID = j(1).layout.sectionId;
x = zeros(numel(j),1);
y = zeros(numel(j),1);
for jix = 1:numel(j)
    jt1(jix) = tile(j(jix));
    tids{jix} = jt1(jix).renderer_id;
    x(jix) = jt1(jix).tform.T(3,1) + jt1(jix).W/2;
    y(jix) = jt1(jix).tform.T(3,2) + jt1(jix).H/2;
end

if nargout>3
    L1 = Msection;
    L1.sectionID = sectionID;
    L1.tiles = jt1;
    L1.X = x(:);
    L1.Y = y(:);
end
