function [mL, tpr, minconf, maxconf] = tile_based_point_pair_errors(mL, A, xout, j, minconf, maxconf)
        % generate point-pair residual information
        res = A*xout;
        resx = zeros(size(res,1)/2,1);
        resy = zeros(size(res,1)/2,1);
        tpr = cell(numel(mL.tiles),1);
        % initialize
        for ix = 1:numel(mL.tiles)
           tpr{ix} = [];
        end
        
        
        % aggregate absolute point-pair error for each tile
        count = 1;
        for ix = 1:size(mL.pm.M,1)
            indx1 = mL.pm.adj(ix,1);
            indx2 = mL.pm.adj(ix,2);
            
            vecx = count:count+mL.pm.np(ix)-1;
            vecy = count+mL.pm.np(ix):count+2*mL.pm.np(ix)-1;
            npoints = numel(vecx);
            rx = sum(abs(res(vecx)),1)/npoints;
            ry = sum(abs(res(vecy)),1)/npoints;
            tpr{indx1} = [tpr{indx1};rx ry];
            tpr{indx2} = [tpr{indx2};rx ry];
            count = count + 2* mL.pm.np(ix);
        end
        
        % average the error
        for ix = 1:numel(mL.tiles)
            mL.tiles(ix).confidence = sum(tpr{ix},1)/size(tpr{ix},1);
            if isempty(tpr{ix})
                mL.tiles(ix).confidence = -50;
            end
        end
        
        
        
        
        %% split into z and display
        ml = split_z(mL);
        [obj, h, rh, A, minconf, maxconf] = show_map_confidence(ml(j), [1], minconf, maxconf);
        
       