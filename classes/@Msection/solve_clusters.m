function [mL, err, R, L_s_mL] = solve_clusters(L_vec, opts, stvec_flag)
% % returns a solved fully assembled coherent group of tiles
% % that could be made of more than one cluster (connected component)
% % 
% % L has full point pair information as obtained for example by using:
% % rc.stack = 'v9_acquire_LC_merged_2';
% % rc.owner='flyTEM';
% % rc.project='FAFB00';
% % rc.baseURL ='http://tem-services.int.janelia.org:8080/render-ws/v1';
% %
% % pm.server = 'http://tem-services.int.janelia.org:8080/render-ws/v1';
% % pm.owner  = 'flyTEM';
% % pm.match_collection = 'v9_2';
% %
% % nfirst = 1502;
% % nlast  = 1502;
% %
% % [L, tIds, PM] = load_point_matches(nfirst, nlast, rc, pm);
% %
% % Author: Khaled Khairy: Janelia Research Campus 2016
% %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3, stvec_flag = 0;end  % in that case we do not assume a starting value, and perform rigid fit

L_s_mL = Msection;
err = [];
R = {};
%parfor_progress(numel(L_vec));
%disp('Solving connected components');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOSI: there is some bug that results in "Conversion to double from Msection is not
% possible." error when the loop is a parfor.
% Also tobd: test accuracy of placement of chuncks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for cix = 1:numel(L_vec)
    %disp(['Processing: ' num2str(cix) ' with #tiles: ' num2str(numel(L_vec(cix).tiles))]);
    if numel(L_vec(cix).tiles)>opts.min_tiles
        
        if stvec_flag==0
            %%% perform rigid transformation
            [ll2r, errR, mL, is, it]  = get_rigid_approximation(L_vec(cix), opts.solver);
            % tiles that are members of reference section (for example that
            % correpsonding to tile 1, and also who did not fail during the
            % solution
            %         indxz = find([ll2r.tiles(:).z]==(ll2r.tiles(1).z) & [ll2r.tiles(:).state]>-2);
            indxz = find([ll2r.tiles(:).z]==(ll2r.tiles(1).z));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % RECOVER POSSIBLE ROTATION AND TRANSLATION RELATIVE TO ORIGINAL CENTER OF MASS
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%% Do rotation first
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Record angle of reference layer before rigid transformation
            anglep = zeros(numel(indxz), 1);
            for ix = 1:numel(indxz)
                anglep(ix) = acosd(L_vec(cix).tiles(indxz(ix)).tform.T(1));
            end
            %%% Record angle of first layer after rigid transformation
            angle = zeros(numel(indxz),1);
            for tix = 1:numel(indxz)
                angle(tix) = acosd(ll2r.tiles(indxz(tix)).tform.T(1));
            end
            % average angle difference for the reference layer
            da = mean(angle-anglep);
            
            %%% calculate center of mass of this connected component
            cm = zeros(numel(ll2r.tiles),2);
            for tix = 1:numel(ll2r.tiles)
                cm(tix,:) = ll2r.tiles(tix).tform.T([3 6]);
            end
            cm = sum(cm)/numel(ll2r.tiles);
            %%% rotate all tiles around this center of mass using the average
            %%% angle da
            for tix = 1:numel(ll2r.tiles)
                ll2r.tiles(tix) = rotate(ll2r.tiles(tix), cm(1), cm(2), da);
            end
            
            %%%%% Now do tanslation for this connected component
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Record tile positions of reference section before rigid transformation
            tp = zeros(numel(indxz),2);%
            for ix = 1:numel(indxz)
                tp(ix,:) = L_vec(cix).tiles(indxz(ix)).tform.T([3 6]);
            end
            %%% Record tile positions of reference section after rigid transformation
            %indxz = find([ll2r.tiles(:).z]==(ll2r.tiles(1).z));
            t = zeros(numel(indxz),2);%
            for ix = 1:numel(indxz)
                t(ix,:) = ll2r.tiles(indxz(ix)).tform.T([3 6]);
            end
            % solve linear system for translation
            b  = tp-t;
            A  = -ones(size(tp,1),1);
            dx = A\b(:,1);        %%recover tile positions (solve for x translation)
            dy = A\b(:,2);        %%recover tile positions (solve for y translation)
            % translate all tiles in the collection by dx and dy
            for tix = 1:numel(ll2r.tiles)
                ll2r.tiles((tix)).tform.T([3 6]) = ll2r.tiles((tix)).tform.T([3 6])  - [dx dy];
            end
        else
            ll2r = L_vec(cix);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if opts.degree==0
            L_s_mL(cix) = ll2r;
        elseif opts.degree==1
            [L_s_mL(cix), err(cix), Res, A, b, B, d, W, K, Lm, xout, LL2, U2, tB, td, invalid] = solve_affine_explicit_region(ll2r, opts);
        elseif opts.degree>1
            [ll3, err(cix), Res, A, b, B, d, W, K, Lm, xout, LL2, U2, tB, td, invalid] = solve_affine_explicit_region(ll2r, opts);
            [L_s_mL(cix), err(cix), Res, A, b, B, d, W, K, Lm, xout, LL2, U2, tB, td] = solve_polynomial_explicit_region(ll3,opts.degree, opts);
        end
        %L_s_mL(cix) = ll;
        R(cix) = {[Res]};
    else
        L_s_mL(cix) = L_vec(cix);
        R(cix) = {[0]};
        err(cix) = nan;
    end
    %parfor_progress;
end
% parfor_progress(0);
%%
mL = concatenate_tiles(L_s_mL, opts.outlier_lambda);