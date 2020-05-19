function [Centres, P] = LineMaxima_v1(Imstack,varargin)
%% centres = LineMaxima_v1(Imstack, varargin)
% Find the centre of the cell in each image of an Imstack by considering
% distance between the row/column maxima in each frame, after applying
% contrast enhancement. 

disp('Finding centre using line maxima method')
tic
NFrames = size(Imstack{1},1);
[NRows, NCols] = size(Imstack{1}{1,1});

P = inputParser();
ParseInputs();

Centres = zeros(2,NFrames);
RowFirsts = zeros(NRows,1);
RowLasts = zeros(NRows,1);
ColFirsts = zeros(NCols,1);
ColLasts = zeros(NCols,1);
NERows = zeros(1,NRows);
NECols = zeros(1,NCols);

if strcmpi(P.Filt,'flat'); Kernel = ones(P.KSize)./P.KSize.^2; end
FiltIm = zeros(NRows,NCols,'like',Imstack{1}{1,1});

for Frame = 1:NFrames    
    ApplyFilt(Imstack);
    % Sharpen image with Laplacian filter - basically a highly tunable
    % way of exaggerating the edges. Time consuming: for an HL60 dataset,
    % 55 s to apply filter to every frame. (55 ms/frame)
    SharpIm = locallapfilt(FiltIm, P.LapSigma, P.LapAlpha, P.LapBeta);
    Thresh = prctile(SharpIm, P.PercentThresh, [1,2]);
    BWIm = SharpIm >= Thresh;
    
    FindFirstsAndLasts()
    Centres(1,Frame) = median(mean([RowFirsts(NERows(NERows~=0)),RowLasts(NERows(NERows~=0))],2));
    Centres(2,Frame) = median(mean([ColFirsts(NECols(NECols~=0)),ColLasts(NECols(NECols~=0))],2));

    MyProgressBar(Frame/NFrames);
end

fprintf('Finished finding centre in %g s\n',toc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ParseInputs
    function ParseInputs()
        FStack = dbstack;
        FName = FStack(2).name;
        fprintf('Parsing inputs to %s\n',FName)
        
        addRequired(P,'Imstack',@(x) ImstackTest(x));
        addParameter(P,'Filt','gaussian',@(x) strcmpi(x,...
            validatestring(x,{'flat','gaussian','binning'},FName,'Filt')))
        addParameter(P,'KSize',4,@(x) validateattributes(x,...
            {'numeric'},{'positive','<',NRows/2,'<',NCols/2},FName,'KSize'))
        addParameter(P,'PercentThresh',95,@(x) validateattributes(x,...
            {'numeric'},{'positive','<',100},FName,'PercentThresh'))
        addParameter(P,'LapSigma',0.1,@(x) validateattributes(x,...
            {'numeric'},{'nonnegative'},FName,'LapSigma'))
        addParameter(P,'LapAlpha',5,@(x) validateattributes(x,...
            {'numeric'},{'positive'},FName,'LapAlpha'))
        addParameter(P,'LapBeta',10,@(x) validateattributes(x,...
            {'numeric'},{'positive'},FName,'LapBeta'))
        
        
        parse(P,Imstack,varargin{:})
        P = P.Results;
        function ImstackTest(x)
            try
                x{1}; %#ok<VUNUS>
            catch
                validateattributes(x,{'cell'},{'nonempty'}, FName, 'Imstack')
            end
            validateattributes(x{1},{'cell'},{'nonempty','ncols',2}, FName, 'Imstack')
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FindFirstsAndLasts
    function FindFirstsAndLasts()
        NERows = (1:NRows) .* (sum(BWIm,2)~=0)';
        NECols = (1:NCols) .* (sum(BWIm,1)~=0);
        for Row = NERows(NERows~=0)
            RowFirsts(Row) = find(BWIm(Row,:),1,'first');
            RowLasts(Row) = find(BWIm(Row,:),1,'last');
        end
        for Col = NECols(NECols~=0)
            ColFirsts(Col) = find(BWIm(:,Col),1,'first');
            ColLasts(Col) = find(BWIm(:,Col),1,'last');
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ApplyFilt
    function ApplyFilt(Imstack)
        % Flatten with ones kernel - maintains an amount of sharpness
        % differently to a Gaussian. I don't know if this is necessary, or
        % if it's better than a Gaussian. You're welcome to try other
        % things
        switch lower(P.Filt)
            case 'flat'
                % Flatten the frame using kernel above
                FiltIm = conv2(Imstack{1}{Frame,1}, Kernel, 'same');
            case 'gaussian'
                % Flatten with gaussian filter
                FiltIm = imgaussfilt(Imstack{1}{Frame,1},(P.KSize-1)/4);
            case 'binning'
                error('Binning not coded yet')
            otherwise
                error('Dafuq has happened in ApplyFilt????')
        end
    end
end