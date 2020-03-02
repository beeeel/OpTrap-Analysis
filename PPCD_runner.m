%% Collection of settings that work for various datasets
function PPCD_runner(varargin)


p = inputParser;

% Which datasets
CellDefault = 'HL60';
DSetsDefault = {'normoxia','drugs','normoxia','hypoxia'} ;
NumsDefault = 1:20;

% What settings
FindVerDefault = 0;
FindDefault = {'Sensitivity',0.95,'Rs',[80 120],'Gfilt',3};
SegDefault = {'iterations', 300, 'method', 'edge','Lsigma',0.1,'Lalpha',5,'Lbeta',10};
UnwrapDefault = {'centering', 1};

% What to save
SummFigDefault = true;
ToSaveDefault = true;
KeepDefault = {'uTaylorParameter','TaylorParameter','Area','centres','filepath'};
SaveAllDefault = false;
            
% Where to save
InfosDirDefault = '/home/ppxwh2/Documents/data/OpTrap/infos/';
FigSaveDirDefault = '~/Documents/data/OpTrap/processing_plots/';

ParseInputs();

for idx = 1:length(p.Results.DSets)
    DSet = p.Results.DSets{idx};
    if idx == 2; p.Results.CellType = 'LS174T'; end
    FolderName = ['/home/ppxwh2/Documents/data/OpTrap/2017_10_movies-from-aishah/' p.Results.CellType '/'];
    
    for Num = p.Results.Nums
        RunNo = num2str(Num);
        if strcmp(p.Results.CellType, 'HL60')
            if strcmp(DSet,'normoxia')
                SetName = 'HL60_normoxia';
                if str2double(RunNo) <= 10
                    FileName = ['100717_' RunNo '_HL60_1.avi'];
                else
                    FileName = ['HL60_' RunNo '_0.020mms-1_1.avi'];
                end
            elseif strcmp(DSet,'drugs')
                SetName = 'HL60_with_drugs';
                FileName = ['190717_HL60_' RunNo '_0.020mm-1_1.avi'];
            end
            Imstack = avi_to_imstack([FolderName SetName '/' FileName]);
        elseif strcmp(p.Results.CellType,'LS174T')
            SetName = [p.Results.CellType '_' DSet];
            FileName = ['200717_' RunNo '_' p.Results.CellType '_' DSet '_1.avi'];
            if strcmp(DSet,'hypoxia')
                FileName(2) = '1';
            end
            Imstack = avi_to_imstack([FolderName FileName]);
        elseif strcmp(p.Results.CellType,'HeLa')
            
        end
        %%
        
        [info, meta] = PostProcessCellDeform_v2(Imstack,'find_cell_v',FindVerDefault,...
            'find_cell',p.Results.FindOpts, 'seg_cell_v',5,'segment_cell', p.Results.SegOpts, ...
            'unwrap_cell_v',2,'unwrap_cell', p.Results.UnwrapOpts);
        
        
        if p.Results.SummFig == true
            MakeSummFig();
        end
        
        if p.Results.ToSave
            if ~p.Results.SaveAll
                FieldNames = fieldnames(info);
                for fld = FieldNames'
                    if sum(strcmp(fld{:},p.Results.KeepFields)) == 0
                        info = rmfield(info, fld{:});
                    end
                end
                save([p.Results.InfosDir 'info_reduced_seg_' SetName '_' FileName(1:end-4) '.mat'], 'info', 'meta');
            else
                save([p.Results.InfosDir 'info_seg_' SetName '_' FileName(1:end-4) '.mat'], 'info', 'meta');
            end
        end
    end
end

    function ParseInputs()
        FName = 'PPCD_runner input validation';
        addParameter(p,'FindVer',FindVerDefault,@(x)validateattributes(x,...
            {'numeric'},{'nonempty','nonnegative','<',3},FName,'FindVer'))
        addParameter(p,'FindOpts',FindDefault,@(x)validateattributes(x,...
            {'cell'},{'nonempty'},FName,'FindOpts'))
        addParameter(p,'SegOpts',SegDefault,@(x)validateattributes(x,...
            {'cell'},{p,'nonempty'},FName,'SegOpts'))
        addParameter(p,'UnwrapOpts',UnwrapDefault,@(x)validateattributes(x,...
            {'cell'},{'nonempty'},FName,'UnwrapOpts'))
        addParameter(p,'KeepFields',KeepDefault,@(x)validateattributes(x,...
            {'cell'},{'nonempty'},FName,'KeepFields'))
        addParameter(p,'ToSave',ToSaveDefault,@(x)validateattributes(x,{'logical'}, {'nonempty'},FName,'ToSave'))
        addParameter(p,'SaveAll',SaveAllDefault,@(x)validateattributes(x,{'logical'},{'nonempty'},FName,'SaveAll'))
        addParameter(p,'SummFig',SummFigDefault,@(x)validateattributes(x,{'logical'},{'nonempty'},FName,'SummFig'))
        addParameter(p,'CellType',CellDefault,@(x)validateattributes(x,...
            {'string','char'},{'nonempty','row','scalartext'},FName,'CellType'))
        addParameter(p,'InfosDir',InfosDirDefault,@(x)validateattributes(x,...
            {'string','char'},{'nonempty','row','scalartext'},FName,'InfosDir'))
        addParameter(p,'FigSaveDir',FigSaveDirDefault,@(x)validateattributes(x,...
            {'string','char'},{'nonempty','row','scalartext'},FName,'FigSaveDir'))
        addParameter(p,'DSets',DSetsDefault,@(x)validateattributes(x,...
            {'cell'},{'nonempty','row'},FName,'DSets'))
        addParameter(p,'Nums',NumsDefault,@(x)validateattributes(x,...
            {'numeric'},{'nonempty','row','nonnegative'},FName,'Nums'))
        
        parse(p,varargin{:});
    end
%%
    function MakeSummFig()
        figure
        subplot(221)
        plot([info.MajorAxisLength]), hold on
        plot([info.MinorAxisLength]), hold off
        legend('Major','Minor')
        title('Regionprops axes')
        subplot(222)
        plot([info.TaylorParameter]), hold on
        plot([info.uTaylorParameter]), hold off
        legend('Regionprops','Unwrap')
        title('Taylor parameters')
        subplot(223)
        plot([info.Area])
        title('Area')
        subplot(224)
        if p.Results.FindVer
            plot([info.radius])
            title('Radius (find_cell)')
        else
            plot([info.mCentres]')
            title('Centres (Line_Maxima)')
        end
        hgsave([p.Results.FigSaveDir, strjoin({SetName,RunNo,'summaryplot'},'_')])
        close

end
end
