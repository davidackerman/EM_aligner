%% Example script to detect level of drift, jumps and warps in collections
%  in comparision to a reference collection.

clc;clear all;kk_clock;

nfirst = 1;
nlast  = 1000;   % preferably set nlast-nfirst = number of CPUs, a multiple of deta_z (below)
delta = 300;  % box size in pixels for x y
delta_z = 16;  % number of z sections spanned by each box --- if using parallel cluster make this number of CPUs
scale = 0.5;  % box will be scaled down by this amount for faster processing
nprobes = 6;   % number of probes (boxes) for each slab

clear rc_vec;

% reference stack: v13 

rc.stack          = ['v13_align'];
rc.owner          ='flyTEM';
rc.project        = 'FAFB00';
rc.service_host   = '10.40.3.162:8080';
rc.baseURL        = ['http://' rc.service_host '/render-ws/v1'];
rc.verbose        = 1;

% test stack v14_align_tps
rc0.stack          = ['v14_align_tps'];
rc0.owner          ='flyTEM';
rc0.project        = 'FAFB00';
rc0.service_host   = '10.40.3.162:8080';
rc0.baseURL        = ['http://' rc0.service_host '/render-ws/v1'];
rc0.verbose        = 1;
rc_vec(1) = rc0;

% test stack v14_align_tps_03
rc1.stack          = ['v14_align_tps_03'];
rc1.owner          ='flyTEM';
rc1.project        = 'FAFB00';
rc1.service_host   = '10.40.3.162:8080';
rc1.baseURL        = ['http://' rc1.service_host '/render-ws/v1'];
rc1.verbose        = 1;
rc_vec(2) = rc1;

%%
[D, I, x1, y1]= detect_jumps_drift_warps(...
    rc, rc_vec, nfirst, nlast, nprobes, delta, delta_z, scale);
% D is a cell array of dimensions nslabs x 1
% nslabs is determined as nfirst:delta_z:nlast
% Also:
% D{1} is of dimensions  (collection number (including reference)) x nprobes
% each element of D{1} is of dimension delta_z x 1 and includes cross
% correlations within that probe
nslabs = numel(D);
dz = 1:delta_z:(nlast-nfirst);dz = dz(end);
%% report correlation results: aggregate
close all;
mean_corr = [];
for six = 1:nslabs
    for pix = 1:nprobes
%         disp('slab  probe collection  z      corr');  % report at most detailed section level
        disp('slab  probe collection  corr');    % report at mean corr per probe level
        for cix = 1:numel(rc_vec)+1
            if cix==1, stack = rc.stack;
            else stack = rc_vec(cix-1).stack;
            end
            im = I{six, cix}{pix}; % image volume for slab six, collection cid and probel pix
            mean_corr(cix, pix) = mean(D{six}{cix, pix});
            str = sprintf('%d\t%d\t%d\t%.4f', six, pix, cix, mean_corr(cix, pix)  );
            disp(str);
        end
        figure; bar(mean_corr);   title(str); drawnow;
    end
end
disp(mean_corr);


%% uncomment to look at the really bad spots
% corr_thresh = 0.15;
% close all;
% for six = 1:nslabs
%     for pix = 1:nprobes
%          disp('slab  probe collection  z      corr');  % report at most detailed section level
%         for cix = 1:numel(rc_vec)+1
%             if cix==1, stack = rc.stack;
%             else stack = rc_vec(cix-1).stack;
%             end
%             im = I{six, cix}{pix}; % image volume for slab six, collection cid and probel pix
%  
%             for zix = 1:dz
%                 %disp([six pix cix zix D{six}{cix, pix}(zix)]);
%                 str = sprintf('%d\t%d\t%d\t%d\t%.4f', six, pix, cix, nfirst+zix-1,...
%                     D{six}{cix, pix}(zix));
%                 disp(str);
%                 if isnan(D{six}{cix, pix}(zix)) || D{six}{cix, pix}(zix)<corr_thresh
%                 figure;imshowpair(im(:,:,zix), im(:,:,zix+1), 'montage');
%                 title(str);
%                 %title([stack ' -- Slab: ' num2str(six) ' probe: ' num2str(pix) ' -- section: ' num2str(nfirst+zix-1) ' and ' num2str(nfirst+zix)]);
%                 pause(1);
%                 end
%             end
%             disp('-----------------------------------');
%         end
%     end
% end




