function [zu, sID, sectionId, z, ns] = get_section_ids(rc, nfirst, nlast)
urlChar = sprintf('%s/owner/%s/project/%s/stack/%s/sectionData', ...
    rc.baseURL, rc.owner, rc.project, rc.stack);
js = webread(urlChar);
sectionId = {js(:).sectionId};
[z, ia]   = sort(([js(:).z]));
sectionId = sectionId(ia);

indx = find(z>=nfirst & z<=nlast);
sectionId = sectionId(indx);% determine the sectionId list we will work with
z         = z(indx);        % determine the zvalues (this is also the spatial order)

% we need unique values, and we need to know how many sectionId's correspond to each unique z value
% usually it is one, but sometimes we have hi/lo dose or other regions
[zu, ia, ic] = unique(z);
count = 1;
sID = {};
for zix = 1:numel(zu)
    ns(zix) =  numel(find(ic==zix));
    vec = {};
    for six = 1:ns(zix)
        vec{six} = sectionId{count};
        sID{zix} = vec;
        count = count + 1;
    end
end