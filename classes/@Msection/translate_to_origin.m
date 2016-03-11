function mL = translate_to_origin(mL)
% mL = get_bounding_box(mL);

for ix = 1:numel(mL.tiles)
    if isa(mL.tiles(ix).tform,  'images.geotrans.PolynomialTransformation2D')
        X(ix) = mL.tiles(ix).tform.A(1);
        Y(ix) = mL.tiles(ix).tform.B(1);
        
    else
        X(ix) = mL.tiles(ix).tform.T(3);
        Y(ix) = mL.tiles(ix).tform.T(6);
        
    end
end

delta = 0;
dx = min(X(:)) + delta;%mL.box(1);
dy = min(Y(:)) + delta;%mL.box(2);
for ix = 1:numel(mL.tiles)
    mL.tiles(ix) = translate_tile(mL.tiles(ix), [dx dy]);
end