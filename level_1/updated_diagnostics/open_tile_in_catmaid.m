function open_tile_in_catmaid(~,~,url,tix)
    if nargin==4
        fprintf('tile index: %d', tix);
    end
    command = sprintf('ssh -X c11u25 "firefox %s & exit"',url);
    system(command);
end