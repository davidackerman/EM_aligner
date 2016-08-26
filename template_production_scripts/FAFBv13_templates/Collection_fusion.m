% test collection fusion function

clc;%clear all;

Slab_definition;
kk_clock;
%%
for ix = 16:nslabs
    

    if ix<2
        error('ix must be larger than or equal to 2');
    elseif ix==2
        %  use the first collection as fixed
        rcfixed                = fine_collection{1};
    else
        rcfixed                = rcfixed_o;
    end
    % configure collection
    rcmoving                = fine_collection{ix};
    % configure output collection
    rcout                  = rcfixed_o;
    overlap = overlapvec(ix-1,:);
    
    
    %%%% just to check
    disp('-----------------------');
    disp('Fixed:');disp(rcfixed);
    disp('Moving:');disp(rcmoving);
    disp('Out:');disp(rcout);
    disp('Overlap:');disp(overlap);
    disp('-----------------------');

    
    %%
    if ix ==2, 
        collection_start = 1;
    else
        collection_start = 0;
    end
    resp = fuse_collections(rcsource, rcfixed, rcmoving, overlap, rcout, collection_start);
    
end
kk_clock;
% delete_renderer_stack(rcout);
