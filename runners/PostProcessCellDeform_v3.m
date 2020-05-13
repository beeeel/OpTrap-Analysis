function [info, meta] = PostProcessCellDeform_v3(Imstack, varargin)
%% [info, meta] = PostProcessCellDeform_v2(Imstack, varargin)
%  Processing stack for OME/TIFF image sets for deformation. 
% Input a cell array from bfopen or avi_to_imstack, and input arguments for
% modules in name-value pairs, e.g.:
%   PostProcessCellDeform_v2(Imstack, 'find_cell', {'Rs',[100, 200]})
% Outputs info struct array with measurements from each frame and meta
% scalar struct with metadata and options
% 
% Branched from v2 during coronavirus times

%% Preamble
% Timing startpoint
StartTime = tic;
% Used by functions
pct = [' ' repmat('%',1,25)];
% Where the file was
FilePath = GetFilePath(Imstack);
% Number of frames
N_frames = size(Imstack{1},1);
% Size of each frame
sz_frame = size(Imstack{1}{1,1});

%% Main body
% Does what it says

PPCD_Par = ParseInputs(Imstack, varargin{:});

[info, meta] = PreallocateInfoMeta(PPCD_Par);

RunFindCell();

RunLineMaxima();

RunUnwrapCell();

RunSegmentCell();

meta.TotalRunTime = toc(StartTime);
%% Get filepath from Imstack
    function FilePath = GetFilePath(Imstack)
        SplitPath = strsplit(Imstack{1}{1,2});
        if length(SplitPath) <= 2
            FilePath = SplitPath{1};
        elseif length(SplitPath) >= 3
            FilePath = strjoin(SplitPath(1:end-1),' '); % Some filenames contain spaces
            warning('Filepath may contain spaces - be careful or my code will give you weird bugs!')
        end
    end
%% Parse input arguments
    function PPCD_Par = ParseInputs(Imstack, varargin)
        % Max version numbers of modules
        ModNames = {'segment_cell','find_cell','unwrap_cell','line_maxima'};
        VMaxFC = 5;
        VMaxSC = 6;
        VMaxUC = 4;
        VMaxLM = 2;
        
        FStack = dbstack;
        FName = FStack(2).name;
        fprintf('Parsing inputs to %s\n',FName)
       
        P = inputParser();
        
        addRequired(P,'Imstack',@(x) ImstackTest(x));
        
        addParameter(P,'segment_cell_v',0, @(x)validateattributes(x,{'numeric'},...
            {'nonnegative','integer','<=',VMaxSC},FName,'segment_cell_v'))
        addParameter(P,'segment_cell_args',{},@(x) ArgsTest(x,ModNames{1}))
        
        addParameter(P,'find_cell_v',0, @(x)validateattributes(x,{'numeric'},...
            {'nonnegative','integer','<=',VMaxFC},FName,'find_cell_v'))
        addParameter(P,'find_cell_args',{},@(x) ArgsTest(x,ModNames{2}))
        
        addParameter(P,'unwrap_cell_v',0, @(x)validateattributes(x,{'numeric'},...
            {'nonnegative','integer','<=',VMaxUC},FName,'unwrap_cell_v'))
        addParameter(P,'unwrap_cell_args',{},@(x) ArgsTest(x,ModNames{3}))
        
        addParameter(P,'line_maxima_v',0, @(x)validateattributes(x,{'numeric'},...
            {'nonnegative','integer','<=',VMaxLM},FName,'line_maxima_v'))
        addParameter(P,'line_maxima_args',{},@(x) ArgsTest(x,ModNames{4}))
        
        parse(P,Imstack,varargin{:})
        PPCD_Par = P.Results;
        
        function ImstackTest(x)
            try
                x{1}; %#ok<VUNUS>
            catch
                validateattributes(x,{'cell'},{'nonempty'}, FName, 'Imstack')
            end
            validateattributes(x{1},{'cell'},{'nonempty','ncols',2}, FName, 'Imstack')
        end
        
        function ArgsTest(x,ArgsName)
            if isempty(x)
                validateattributes(x,{'cell'},{'2d'},FName, ArgsName)
            else
                validateattributes(x,{'cell'},{'row'},FName, ArgsName)
            end
        end
    end

%% Preallocate structs for output
    function [info, meta] = PreallocateInfoMeta(Par)
        % The order of fields in {regionFields, fields} must be the same as
        % the order of defaults in values. Hopefully the linebreaks help?
        RegionFields = {'MinorAxisLength', 'MajorAxisLength', 'Orientation', ...
            'Centroid', 'Perimeter', 'Area', 'Eccentricity'};
        RegionValues = {0,        0,      0, ...
            zeros(1,2),         0,      0,      0,};
        
        % Concatenating cell arrays with [] is quicker than deconstructing with {:}
        % and reconstructing inside {...}
        SegFields = [RegionFields, 'mask', 'seg_fails','ellipse_fits', 'TaylorParameter', 'Flatness'];
        SegValues = [RegionValues, false(sz_frame), uint8(0), zeros(6,1), 0, 0];
        
        FindFields = {'centres', 'radius'};
        FindValues = {[0;0], 0};
        
        UnwrapFields = {'uMajorAxisLength','uMinorAxisLength',...
            'uOrientation','uTaylorParameter','uFlatness','uOffset','uFitErrs'};
        UnwrapValues = {0, 0, ...
            0, 0, 0, zeros(3,1), zeros(3,1)};
        
        LineMaxFields = {'mCentres', 'mWidths'};
        LineMaxDefaults = {[0;0], [0;0]};
        
        InfoFields = [FindFields, UnwrapFields, SegFields, LineMaxFields, 'filepath'];
        InfoValues = [FindValues, UnwrapValues, SegValues, LineMaxDefaults, FilePath{1}];
        
        MetaFields = {'filepath', 'N_Frames','Frame_size', 'TotalRunTime', ...
            'find_cell_time', 'unwrap_cell_time', 'segment_cell_time', 'line_maxima_time',...
            'find_cell_args', 'find_cell_v', 'segment_cell_args', 'segment_cell_v', ...
            'unwrap_cell_args', 'unwrap_cell_v', 'line_maxima_args', 'line_maxima_v'};
        MetaValues = {FilePath{1}, N_frames, sz_frame,0, ...
            0, 0, 0, 0,...
            Par.find_cell_args, Par.find_cell_v, Par.segment_cell_args, Par.segment_cell_v, ...
            Par.unwrap_cell_args, Par.unwrap_cell_v, Par.line_maxima_args, Par.line_maxima_v};
        
        % Create one struct for info and one for metadata
        info(1:N_frames) = cell2struct(InfoValues, InfoFields, 2);
        meta = cell2struct(MetaValues, MetaFields, 2);
        
        % Tidy it to remove excess fields
        if Par.segment_cell_v == 0
            info = rmfield(info,SegFields );
            meta = rmfield(meta, {'segment_cell_time','segment_cell_args'});
        end
        if Par.find_cell_v == 0
            info = rmfield(info, FindFields);
            meta = rmfield(meta, {'find_cell_time','find_cell_args'});
        end
        if Par.unwrap_cell_v == 0
            info = rmfield(info, UnwrapFields);
            meta = rmfield(meta, {'unwrap_cell_time','unwrap_cell_args'});
        end
        if Par.line_maxima_v == 0
            info = rmfield(info, LineMaxFields);
            meta = rmfield(meta, {'line_maxima_time','line_maxima_args'});
        end
    end

%% Run find_cell
    function RunFindCell()
        if PPCD_Par.find_cell_v ~= 0
            before = toc(StartTime);
            % 1-Find cells,
            fprintf('%s\nFinding cells\nUsing find_cell version %i\n%s\n',...
                pct, PPCD_Par.find_cell_v, pct)
            
            % Starting cell finding - returns a struct array with info on the best
            % circle in each image, tracking of when it failed, and a crop array
            % showing where to crop for the best circle in each image.
            switch PPCD_Par.find_cell_v
                case 1
                    cell_dat = find_cell(Imstack, PPCD_Par.find_cell_args{:});
                case 2
                    [cell_dat, FCPar] = find_cell_v2(Imstack, PPCD_Par.find_cell_args{:});
                case 5
                    [Centres, Widths, FCPar] = find_cell_v5(Imstack, PPCD_Par.find_cell_args{:});
            end
            if PPCD_Par.find_cell_v <= 3
                % Transfer the info from cell_dat struct into info struct
                InfoFields = fieldnames(cell_dat);
                for frame = 1:N_frames
                    for f_no = 1:numel(InfoFields)
                        info(frame).(InfoFields{f_no}) = cell_dat(frame).(InfoFields{f_no});
                    end
                end
            else
                for frame = 1:N_frames
                    info(frame).centres = Centres(:,frame);
                    info(frame).radius = Widths(:,frame);
                end
            end
            % Save parameters from find_cell
            meta.find_cell_args = FCPar;
            
            clear cell_dat FCPar Centres Widths
            meta.find_cell_time = toc(StartTime) - before;
            fprintf('%s\nFound cells\n%g s elapsed\n%s\n',...
                pct, meta.find_cell_time, pct)
        else
            fprintf('%s\nSkipping find_cell\n%s\n', pct, pct)
        end
    end

%% Run LineMaxima
    function RunLineMaxima()
        if PPCD_Par.line_maxima_v ~= 0
            fprintf('%s\nFinding cells\nUsing line_maxima version %i\n%s\n',...
                pct, PPCD_Par.line_maxima_v, pct)
            before = toc(StartTime);
            switch PPCD_Par.line_maxima_v
                case 1
                    [Centres, LMPar] = LineMaxima_v1(Imstack,PPCD_Par.line_maxima_args{:});
                case 2
                    [Centres, Widths, LMPar] = LineMaxima_v2(Imstack,PPCD_Par.line_maxima_args{:});
                otherwise
                    error('huhnknown line_maxima_v')
            end
            % Parse outputs
            for frame = 1:N_frames
                info(frame).mCentres = Centres(:,frame);
                info(frame).mWidths = Widths(:,frame);
            end
            meta.line_maxima_args = LMPar;
            meta.line_maxima_time = toc(StartTime) - before;
            fprintf('%s\nFound cells\n%g s elapsed\n%s\n',...
                pct, meta.line_maxima_time, pct)
            clear Centres LMPar
        else
            fprintf('%s\nSkipping line_maxima\n %s\n', pct, pct)
        end
    end

%% Run unwrap_cell
    function RunUnwrapCell()
        if PPCD_Par.unwrap_cell_v ~= 0
            fprintf('%s\nUnwrapping cells\nUsing unwrap_cell version %i\n%s\n',...
                pct, PPCD_Par.unwrap_cell_v, pct)
            before = toc(StartTime);
            
            % Take cell location dependent on how it was found
            if PPCD_Par.find_cell_v == 0 && PPCD_Par.line_maxima_v == 0
                error('You need to set a version for find_cell or line_maxima in order to run unwrap_cell')
            elseif PPCD_Par.find_cell_v ~= 0
                Centres = [info.centres];
                Radii = [info.radius];
            else
                Centres = [info.mCentres];
                Radii = repmat(min(sz_frame)/2,1,N_frames);
            end
            
            % Do the unwrapping
            switch PPCD_Par.unwrap_cell_v
                case 1
                    UnwrapFits = unwrap_cell_v1(Imstack, Centres, Radii, ...
                        PPCD_Par.unwrap_cell_args{:});
                case 2
                    [UnwrapFits, ~, ~, ~, UnwrapOffset, FitErrs, UCPar] = ...
                        unwrap_cell_v2(Imstack, Centres, Radii, PPCD_Par.unwrap_cell_args{:});
                case 4
                    [UnwrapFits, ~, ~, UnwrapOffset, FitErrs, UCPar] = unwrap_cell_v4(Imstack, Centres, Radii, PPCD_Par.unwrap_cell_args{:});
                otherwise
                    error('huh')
            end
            
            % Parse outputs
            for frame = 1:N_frames
                info(frame).uMajorAxisLength = UnwrapFits(1,frame);
                info(frame).uMinorAxisLength = UnwrapFits(2,frame);
                info(frame).uOrientation = UnwrapFits(3,frame);
                info(frame).uFitErrs = FitErrs(1:3,frame);
                
                % Taylor's deformation parameter : (LongestAxisLength-ShortestAxisLength)/(LongestAxisLength+ShortestAxisLength)
                info(frame).uTaylorParameter = ( UnwrapFits(1,frame) - ...
                    UnwrapFits(2,frame)) / ( UnwrapFits(1,frame) ...
                    + UnwrapFits(2,frame));
                
                % Flatness : (majorAxisLength-MinorAxisLength)/MajorAxisLength
                info(frame).uFlatness = ( UnwrapFits(1,frame) - ...
                    UnwrapFits(2,frame)) /  UnwrapFits(1,frame);
            end
            % Offset is not always found - depends on args. This hotfix is
            % ugly af but it works
            if ~isempty(UnwrapOffset)
                for frame = 1:N_frames
                    info(frame).uOffset = UnwrapOffset(:,frame);
                end
            end
            meta.unwrap_cell_args = UCPar;
            meta.unwrap_cell_time = toc(StartTime) - before;
            
            fprintf('%s\nUnwrapped cells\n%g s elapsed\n%s\n',...
                pct, meta.unwrap_cell_time, pct)
            clear Centres Radii UnwrapFits UnwrapOffset FitErrs UCPar
        else
            fprintf('%s\nSkipping unwrap_cell\n %s\n', pct, pct)
        end
    end
%% Run segment_cell
    function RunSegmentCell()
        if PPCD_Par.segment_cell_v ~= 0
            fprintf('%s\nMasking cells\nUsing segment_cell version %i\n%s\n',...
                pct, PPCD_Par.segment_cell_v, pct)
            before = toc(StartTime);
            %% Run segment_cell
            % Segment cell to get masks array - segment_cell_v2 and newer iterate over
            % image stack.
            switch PPCD_Par.segment_cell_v
                case 1
                    % V1: Use edge detection, dilation and erosion
                    % I've not had success using this
                    masks = zeros([sz_frame, N_frames],'uint8');
                    for frame = 1:N_frames
                        masks(:,:,frame) = segment_cell(Imstack{1}{frame,1});
                    end
                case 2
                    % V2: Use Canny edge detection, followed by thresholding and filtering
                    % This only works for very nice looking cells
                    masks = segment_cell_v2(Imstack, 'crop', [info.crop],  ...
                        'fails', [info.find_fails], PPCD_Par.segment_cell_args{:});
                case 3
                    % V3: Use Laplacian filtering for edge enhancement, thresholding and
                    % dilation.
                    % This is more robust than v2.
                    masks = segment_cell_v3(Imstack, 'crop', [info.crop], PPCD_Par.segment_cell_args{:});
                case 4
                    masks = segment_cell_v4(Imstack);
                case 5
                    if PPCD_Par.find_cell_v ~=0
                        [masks, fits, SegFails] = segment_cell_v5(Imstack, 'crop', [info.crop], ...
                            'radius', [info.radius], 'fails', [info.find_fails],...
                            'centres', [info.centres], PPCD_Par.segment_cell_args{:});
                    else
                        [masks, fits, SegFails, SCPar] = segment_cell_v5(Imstack, PPCD_Par.segment_cell_args{:});
                    end
                case 6
                    [masks, fits, SegFails] = segment_cell_v6(Imstack, PPCD_Par.segment_cell_args{:});
                otherwise
                    error('Unknown segment_cell version %d', segment_cell_v);
            end
            
            
            %% Put masks into info struct
            for frame = 1:N_frames
                info(frame).mask = masks(:,:,frame);
                info(frame).seg_fails = SegFails(frame);
                if PPCD_Par.segment_cell_v == 5
                    info(frame).ellipse_fits = fits(frame, :)';
                end
            end
            meta.segment_cell_args = SCPar;
            fprintf('%s\nMasked cells\n%g s elapsed\n%s\nMeasuring Cells\n%s\n',...
                pct, toc(StartTime)-before, pct, pct)
            %% Get region properties
            for frame = 1 : N_frames
                % (optional:) Only look at the cells that were found - if you change
                % the number below to 1, it will skip frames where find_cell failed.
                if 1% info(frame).find_fails ~= -1
                    % Use Regionprops to extract information and populate info structure
                    props = regionprops(imfill(info(frame).mask,'holes'),RegionFields);
                    % Put the properties into the info struct - iterate over the list
                    % of field names.
                    % If the mask is empty, props is a 0x1 struct array, which creates
                    % errors for dot indexing. In this case, put NaNs into info.
                    
                    % If regionprops finds multiple regions in the image (e.g.: other
                    % cells in FoV), the distance to the centre for each centroid is
                    % calculated. The object closest to the centre is chosen.
                    Names = fieldnames(props);
                    if size(props,1) > 0
                        % Store the number of found objects
                        info(frame).NRegions = size(props,1);
                        % If only one object is found, just take it
                        if size(props,1) == 1
                            for field = 1:numel(Names)
                                info(frame).(Names{field}) = props.(Names{field});
                            end
                        else
                            % If multiple are found, calculate the distances
                            dists = zeros(size(props));
                            for idx = 1:size(props,1)
                                dists(idx) = sum((props(idx).Centroid - size(Imstack{1}{1,1}/2)).^2);
                            end
                            % Choose the closest to the centre and take it
                            [~, idx] = min(dists);
                            for field = 1:numel(Names)
                                info(frame).(Names{field}) = props(idx).(Names{field});
                            end
                        end
                    else
                        % If none are found, store NaNs
                        for field = 1:numel(Names)
                            if ~strcmp(Names{field},'Centroid')
                                info(frame).(Names{field}) = NaN;
                            else
                                info(frame).(Names{field}) = [NaN, NaN];
                            end
                        end
                    end
                    
                    % Taylor's deformation parameter : (LongestAxisLength-ShortestAxisLength)/(LongestAxisLength+ShortestAxisLength)
                    info(frame).TaylorParameter = (info(frame).MajorAxisLength - ...
                        info(frame).MinorAxisLength) / (info(frame).MajorAxisLength ...
                        + info(frame).MinorAxisLength);
                    
                    % Flatness : (majorAxisLength-MinorAxisLength)/MajorAxisLength
                    info(frame).Flatness = (info(frame).MajorAxisLength - ...
                        info(frame).MinorAxisLength) / (info(frame).MajorAxisLength);
                    
                else
                    % When no circle is detected, there's no mask, and so nothing we
                    % can do :(
                    info(frame).MinorAxisLength=nan;
                    info(frame).MajorAxisLength=nan;
                    info(frame).Orientation=nan;
                    info(frame).centroid=nan;
                    info(frame).Perimeter=nan;
                    info(frame).Area=nan;
                    info(frame).Eccentricity=nan;
                    info(frame).TaylorParameter=nan;
                    info(frame).Flatness=nan;
                end
                prog = ceil(100 * frame / N_frames);
                fprintf('%s\r',['[' repmat('=',1,prog) repmat(' ',1,100-prog) ']'])
                
            end
            fprintf('%s\r',repmat(' ',1,104))
            clear masks Imstack fits SCPar SegFails
            meta.segment_cell_time = toc(StartTime) - before;
            fprintf('%s\nSegmented cells\n%g s elapsed\n%s\n',...
                pct, meta.segment_cell_time, pct)
        else
            fprintf('%s\nSkipping segment_cell\n%s\n',...
                pct, pct)
        end
    end
end