function obj = translate_to_origin(obj)
% mL = get_bounding_box(mL);

for ix = 1:numel(obj.tiles)
    if isa(obj.tiles(ix).tform,  'images.geotrans.PolynomialTransformation2D')
        X(ix) = obj.tiles(ix).tform.A(1);
        Y(ix) = obj.tiles(ix).tform.B(1);
        
    else
        X(ix) = obj.tiles(ix).tform.T(3);
        Y(ix) = obj.tiles(ix).tform.T(6);
        
    end
end

delta = 0;
dx = min(X(:)) + delta;%mL.box(1);
dy = min(Y(:)) + delta;%mL.box(2);
for ix = 1:numel(obj.tiles)
    obj.tiles(ix) = translate_tile(obj.tiles(ix), [dx dy]);
end