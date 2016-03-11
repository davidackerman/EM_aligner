classdef tile
    properties
        z;
        sectionId;
        id;
        renderer_id = '0';
        tform = affine2d;      % transformation
        col;
        row;
        cam;
        path;
        mask = [];
        temca_conf = -999;
        rot;
        H = 0;
        W = 0;
        state = 1;             % 0 = this tile is marked for removal and should not be rendered (or has invalid transformation information)
        owner = 'flyTEM';
        project='FAFB00';
        stack = '';
        server = 'http://10.40.3.162:8080/render-ws/v1'; %'http://tem-services.int.janelia.org:8080/render-ws/v1'; % 
        features = [];        % store features information
        validPoints = [];     % point locations corresponding to features
        featuresMethod = 'SURF'; % method to be used for calculating features
        % surf paramters
        SURF_NumOctaves = 2;
        SURF_NumScaleLevels = 3;
        SURF_MetricThreshold = 2000;
        fetch_local = 0;
        
    end
    
    
    methods
        function obj = tile(z, id, t1, t2, t3, t4, t5, t6,...
                col, row, cam, path, temca_conf, rotation, renderer_id)
            if (nargin==1 & isstruct(z))   % we are reading from a struct that was produced from json from Renderer db
                p = z;
                obj.z = p.z;
                obj.id = -999;
                obj.renderer_id = p.tileId;
                obj.H = p.height;
                obj.W = p.width;
                obj.cam = str2double(p.layout.camera);
                obj.col = p.layout.imageCol;
                obj.row = p.layout.imageRow;
                obj.rot = p.layout.rotation;
                obj.temca_conf = str2double(p.layout.temca);
                obj.sectionId = p.layout.sectionId;
                
                if strcmp(p.mipmapLevels.x0.imageUrl(1:4), 'file')
                    obj.path = p.mipmapLevels.x0.imageUrl(6:end);
                end
                if strcmp(p.mipmapLevels.x0.imageUrl(1:4), 'file')
                    obj.mask = p.mipmapLevels.x0.maskUrl(6:end);
                end
                % check if we have a list of transformations, in which case
                % we need
                if numel(p.transforms.specList)>1,
                    p.transforms.specList = p.transforms.specList(end);
                end
                
                Tstr = p.transforms.specList.dataString;
                Tdouble    = str2double(strsplit(Tstr));

                if strcmp(p.transforms.specList.className, 'mpicbg.trakem2.transform.AffineModel2D')
                    T(3,3) = 1;
                    T(1,1) = Tdouble(1);
                    T(2,1) = Tdouble(2);
                    T(2,1) = Tdouble(3);
                    T(2,2) = Tdouble(4);
                    T(3,1) = Tdouble(5);
                    T(3,2) = Tdouble(6);
                    obj.tform.T = T;
                else
                    warning('only affine implemented so far for reading from Renderer');
                end
            end
            if nargin>1
                obj.z = z;
                obj.id = id;
                T = zeros(3,3);
                T(3,3) = 1;
                T(1,1) = t1;
                T(2,1) = t2;
                T(3,1) = t3;
                T(1,2) = t4;
                T(2,2) = t5;
                T(3,2) = t6;
                obj.tform.T = T;
                %                 obj.tform.T(1,1) = t1;
                %                 obj.tform.T(2,1) = t2;
                %                 obj.tform.T(3,1) = t3;
                %                 obj.tform.T(1,2) = t4;
                %                 obj.tform.T(2,2) = t5;
                %                 obj.tform.T(3,2) = t6;
                obj.col = col;
                obj.row = row;
                obj.cam = cam;
                obj.path = path;
                obj.temca_conf = -999;
                obj.rot = 0;
                if nargin>12
                obj.temca_conf = temca_conf;
                end
                if nargin>13
                obj.rot = rotation;
                end
                if nargin>14
                    obj.renderer_id = renderer_id;
                end
                %if ~exist(obj.path,'file'), disp(['Tile not found: ' obj.path]);end
            end
        end
        function im = get_mask(obj)
            try
            im = [];
            if ~isempty(obj.mask)
                im = imread(obj.mask);
            end
            catch err_mask
                disp('Error retrieving mask');
            end
        end
        function im = get_image(obj, filter)
            if obj.fetch_local
                im = imread(obj.path);
            else
              im = get_image_renderer(obj, 1.0, 'true');
            end
%             if strcmp(obj.path(1:4), 'http')
%                 im = get_image_renderer(obj, 1.0, '"true"');
%             elseif ~exist(obj.path,'file'),
%                 disp(['Tile not found: ' obj.path]);
%                 im = [];
%             else  %%% then obj.path is just a path to a file
%                 im = imread(obj.path);
%             end
            if nargin==2 && ~isempty(im)
                if strcmp(filter, 'bkgrd1')
                    im = background_filter(im,1);
                elseif strcmp(filter, 'bkgrd2')
                    im = background_filter(im,2);
                elseif strcmp(filter, 'bkgrd3')
                    im = background_filter(im,3);
                elseif strcmp(filter, 'histeq')
                    im = histeq(im);
                end
            end
            mask = get_mask(obj);
            im(~mask) = 0;
        end
        
        function im = get_image_renderer(obj, scale, filter)
            % filter must be a string, either "true" or "false"
            %url = 'http://tem-services.int.janelia.org:8080/render-ws/v1/owner/flyTEM/project/FAFB00/stack/v5_acquire/tile/150127175351050044.3826.0/scale/1.0/jpeg-image?filter=true';
            url = sprintf('%s/owner/%s/project/%s/stack/%s/tile/%s/jpeg-image?scale=%s&filter=%s',...
                obj.server, obj.owner, obj.project, obj.stack, obj.renderer_id, num2str(scale), filter);
%            url = sprintf('%s/owner/%s/project/%s/stack/%s/tile/%s/jpeg-image?scale=%s', obj.server, obj.owner, obj.project, obj.stack, obj.renderer_id, num2str(scale));

            options = weboptions('Timeout', 60);
            try
               im = webread(url, options);
           catch err_ip_address
               pause(1);
               try
               im = webread(url,options);
               catch err_ip_address2
                   pause(2);
                   im = webread(url,options);
               end
           end
            %im = webread(url, options);
            im = rgb2gray(im);
        end
        
        function [im,R] = get_warped_image(obj)
            %warning('Obsolete: using get_image_renderer');
            R = [];
            if obj.fetch_local
                if ~exist(obj.path,'file'),
                    disp(['Tile not found: ' obj.path]);
                    im = [];
                else
                    
                    %%% deal with case when tile is polynomial
                    if strcmp(class(obj.tform), 'images.geotrans.PolynomialTransformation2D') ...
                            ||strcmp(class(obj.tform), 'nonlin2d')
                        px = obj.W/2;
                        py = obj.H/2;
                        %                     U = transformPointsInverse(obj.tform,[px py]);
                        %                     X = U(1);
                        %                     Y = U(2);
                        
                        bo = [0 0;obj.W 0;0 obj.H;obj.W obj.H];
                        p = transformPointsInverse(obj.tform,bo);
                        Wbox = [min(p(:,1)) min(p(:,2)) max(p(:,1))-min(p(:,1)) max(p(:,2))-min(p(:,2))];
                        disp('Warning: Scale is set at 0.2');
                        [im] =  render_poly_03(obj, 0.2, Wbox, 0, 0, [],[]);
                        %                    figure;imshow(im);
                        
                    else
                        [im,R] = imwarp(imread(obj.path), obj.tform, 'Interp', 'cubic');
                    end
                end
            else
                im = get_image_renderer(obj, 1.0, 'true');
            end
        end
        
        function im = show(obj, filter)
            if nargin<2, filter = 'none';end
            warning off;
            disp('Raw image (no warping) is being rendered');
            im = get_image(obj, filter);
            imshow(im);
            warning on;
        end
        function im = show_warped(obj)
             warning off;
            disp('Warped image (transformed according to tform) is being rendered');
            im = get_warped_image(obj);
            imshow(im);
            warning on;
        end
        function F = show_fft(obj)
            F = fftshift(fft2(get_image(obj))); % Center FFT
            F = abs(F); % Get the magnitude
            F = log(F+1); % Use log, for perceptual scaling, and +1 since log(0) is undefined
            F = mat2gray(F); % Use mat2gray to scale the image between 0 and 1
            imshow(F,[]); % Display the result
        end
        function obj = set_info(obj)
            if obj.H==0,
                if ~isempty(obj.path)
                info = imfinfo(obj.path);     % slow
                obj.H = info.Height;
                obj.W = info.Width;
                end
            end
        end
    end
end