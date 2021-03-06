%% Collection of settings that work for various datasets
function PPCD_runner(varargin)

Par = [];
% Which datasets
CellDefault = {'LS174T','HL60','MV411'};
DSetsDefault = {{'normoxia','hypoxia'},{'normoxia','with_drugs'},{'normoxia','with_drugs'}};
NumsDefault = {{1:20, 1:20},{1:30, 1:20},{1:30, 1:20}};

% What settings
FindVerDefault = 5;
FindDefault = {};
SegVerDefault = 0;
SegDefault = {};% {'iterations', 300, 'method', 'edge','Lsigma',0.1,'Lalpha',5,'Lbeta',10};
UnwrapVerDefault = 0;%4;
UnwrapDefault = {};

% What to save
SummFigDefault = true;
ToSaveDefault = true;
KeepDefault = {'uTaylorParameter','uMinorAxisLength','uMajorAxisLength','uOrientation','uOffset','uFitErrs','centres','radius','filepath'};
SaveAllDefault = false;
            
% Where to save
InfosDirDefault = '~/Documents/data/OpTrap/infos/';
FigSaveDirDefault = '~/Documents/data/OpTrap/processing_plots/';
SuffixDefault = 'find5';

ParseInputs();

for CTidx = 1:length(Par.CellType)
    FolderName = ['~/Documents/data/OpTrap/2017_10_movies-from-aishah/' Par.CellType{CTidx} '/'];
    disp(FolderName)
    for Didx = 1:length(Par.DSets{CTidx})
        DSet = Par.DSets{CTidx}{Didx};
        disp(DSet)
        %%{
        for Num = Par.Nums{CTidx}{Didx}
            NumStr = num2str(Num);
            disp(NumStr)
            Imstack = LoadImstackInfoMeta(Par.CellType{CTidx}, DSet, NumStr, true);
            %%
            [info, meta] = PostProcessCellDeform_v3(Imstack,...
                'find_cell_v',Par.FindVer, 'find_cell_args',Par.FindOpts, ...
                'segment_cell_v',Par.SegVer,'segment_cell_args', Par.SegOpts,...
                'unwrap_cell_v',Par.UnwrapVer,'unwrap_cell_args', Par.UnwrapOpts);
            
            if Par.SummFig == true
                MakeSummFig(info);
            end
            
            if Par.ToSave
                if ~Par.SaveAll
                    FieldNames = fieldnames(info);
                    for fld = FieldNames'
                        if sum(strcmp(fld{:},Par.KeepFields)) == 0
                            info = rmfield(info, fld{:});
                        end
                    end
                    save([Par.InfosDir 'info_reduced_' strjoin({Par.CellType{CTidx}, DSet, NumStr, Par.Suffix},'_') '.mat'], 'info', 'meta');
                else
                    save([Par.InfosDir 'info_' strjoin({Par.CellType{CTidx}, DSet, NumStr, Par.Suffix},'_') '.mat'], 'info', 'meta');
                end
            end
        end
        %}
    end
end

    function ParseInputs()
        p = inputParser();
        FName = 'PPCD_runner input validation';
        addParameter(p,'FindVer',FindVerDefault,@(x)validateattributes(x,...
            {'numeric'},{'nonempty','nonnegative'},FName,'FindVer'))
        addParameter(p,'FindOpts',FindDefault,@(x)validateattributes(x,...
            {'cell'},{'nonempty'},FName,'FindOpts'))
        addParameter(p,'SegVer',SegVerDefault,@(x)validateattributes(x,...
            {'numeric'},{'nonempty','nonnegative'},FName,'SegVer'))
        addParameter(p,'SegOpts',SegDefault,@(x)validateattributes(x,...
            {'cell'},{p,'nonempty'},FName,'SegOpts'))
        addParameter(p,'UnwrapVer',UnwrapVerDefault,@(x)validateattributes(x,...
            {'numeric'},{'nonempty','nonnegative'},FName,'UnwrapVer'))
        addParameter(p,'UnwrapOpts',UnwrapDefault,@(x)validateattributes(x,...
            {'cell'},{'nonempty'},FName,'UnwrapOpts'))
        addParameter(p,'KeepFields',KeepDefault,@(x)validateattributes(x,...
            {'cell'},{'nonempty'},FName,'KeepFields'))
        addParameter(p,'ToSave',ToSaveDefault,@(x)validateattributes(x,{'logical'}, {'nonempty'},FName,'ToSave'))
        addParameter(p,'SaveAll',SaveAllDefault,@(x)validateattributes(x,{'logical'},{'nonempty'},FName,'SaveAll'))
        addParameter(p,'SummFig',SummFigDefault,@(x)validateattributes(x,{'logical'},{'nonempty'},FName,'SummFig'))
        addParameter(p,'CellType',CellDefault,@(x)validateattributes(x,...
            {'cell'},{'nonempty','row'},FName,'CellType'))
        addParameter(p,'InfosDir',InfosDirDefault,@(x)validateattributes(x,...
            {'string','char'},{'nonempty','row','scalartext'},FName,'InfosDir'))
        addParameter(p,'FigSaveDir',FigSaveDirDefault,@(x)validateattributes(x,...
            {'string','char'},{'nonempty','row','scalartext'},FName,'FigSaveDir'))
        addParameter(p,'Suffix',SuffixDefault,@(x)validateattributes(x,...
            {'string','char'},{'nonempty','row','scalartext'},FName,'Suffix'))
        addParameter(p,'DSets',DSetsDefault,@(x)ValidateDSets(x))
        addParameter(p,'Nums',NumsDefault,@(x)validateattributes(x,...
            {'cell'},{'nonempty','row'},FName,'Nums'))
        
        parse(p,varargin{:});
        Par = p.Results;
        function ValidateDSets(x)
            validateattributes(x,{'cell'},{'nonempty','row'},FName,'DSets')
            for y = x
                validateattributes(y,{'cell'},{'nonempty'},FName,'DSets')
                for z = y
                    validateattributes(z,{'string','char'},{'nonempty','scalartext'},FName,'DSets')
                end
            end
        end
    end
%%
    function Imstack = N_LoadImstack()
        %{
            if strcmp(Par.CellType{CTidx}, 'HL60')
                if strcmp(DSet,'normoxia')
                    SetName = 'HL60_normoxia';
                    if str2double(RunNo) <= 10
                        FileName = ['100717_' RunNo '_HL60_1.avi'];
                    else
                        FileName = ['HL60_' RunNo '_0.020mms-1_1.avi'];
                    end
                elseif strcmp(DSet,'with_drugs')
                    SetName = 'HL60_with_drugs';
                    FileName = ['190717_HL60_' RunNo '_0.020mm-1_1.avi'];
                end
                Imstack = avi_to_imstack([FolderName SetName '/' FileName]);
            elseif strcmp(Par.CellType{CTidx},'LS174T')
                SetName = [Par.CellType{CTidx} '_' DSet];
                FileName = ['200717_' RunNo '_' Par.CellType{CTidx} '_' DSet '_1.avi'];
                if strcmp(DSet,'hypoxia')
                    FileName(2) = '1';
                end
                Imstack = avi_to_imstack([FolderName FileName]);
            elseif strcmp(Par.CellType{CTidx},'HeLa')
            elseif strcmp(Par.CellType{CTidx}, 'MV411')
                if strcmp(DSet,'normoxia')
                    SetName = 'MV411_normoxia';
                    if str2double(RunNo) <= 10
                        FileName = ['100717_' RunNo '_ MV411_1.avi'];
                    elseif str2double(RunNo) ~= 25
                        FileName = ['MV411_' RunNo '_0.020mms-1_1.avi'];
                    else
                        FileName = ['MV411_' RunNo '_0.020mms-1_2.avi'];
                    end
                elseif strcmp(DSet,'with_drugs')
                    SetName = 'MV411_with_drugs';
                    FileName = ['180717_' RunNo '_MV411_0.020mms-1_1.avi'];
                end
                Imstack = avi_to_imstack([FolderName SetName '/' FileName]);
                
            end
        %}
    end
%%
    function MakeSummFig(info)
        NX = 3;
        NY = 2;
        figure
        subplot(NX, NY, 1)
        cla, hold on
        if meta.segment_cell_v
            plot([info.MajorAxisLength])
            plot([info.MinorAxisLength])
            title('Regionprops axes')
        elseif meta.unwrap_cell_v
            plot([info.uMajorAxisLength])
            plot([info.uMinorAxisLength])
            title('Unwrap axes')
        end
        legend('Major','Minor')
        subplot(NX, NY, 2)
        cla, hold on
        if meta.segment_cell_v && meta.unwrap_cell_v
            plot([info.TaylorParameter])
            plot([info.uTaylorParameter])
            legend('Regionprops','Unwrap')
        elseif meta.unwrap_cell_v
            plot([info.uTaylorParameter])
            legend('Unwrap')
        elseif meta.segment_cell_v
            plot([info.TaylorParameter])
            legend('Regionprops')
        end
        title('Taylor parameters')
        subplot(NX, NY, 3)
        if meta.segment_cell_v
            plot([info.Area])
        end
        title('Area')
        subplot(NX, NY, 4)
        if Par.FindVer
            plot([info.radius])
            title('Radius (find\_cell)')
        else
            plot([info.mCentres]')
            title('Centres (Line\_Maxima)')
        end
        if NY == 3
            subplot(NX, NY, 5)
            if isfield(info,'mCentres')
                plot([info.mCentres]')
                title('Centres (Line\_Maxima)')
            end
            subplot(NX, NY, 6)
            if isfield(info,'centres')
                plot([info.centres]')
                title('Centres (find\_cell)')
            end
        end
        hgsave([Par.FigSaveDir, strjoin({Par.CellType{CTidx},DSet,NumStr,'summaryplot'},'_')])
        close
        
    end
end
