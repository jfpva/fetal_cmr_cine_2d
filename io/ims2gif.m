function [filepath] = ims2gif(ims,varargin)
%IMS2GIF Create animated gif from MR image sets.
%
%   IMS2GIF(IMS) creates an animated gif from 3D array IMS. Images are 
%   taken from the first two dimensions in IMS with the third dimension 
%   assumed to be time or slice. 
%
%   IMS2GIF(IMS,'Option',OptionValue,...) sets specified options. 
%
%   Options inlcude:
%
%       'filename':         name of .gif
%       'filedir':          directory to save .gif
%       'imlimits':         image limits, [min, max]
%       't':                vector of times for each image in milliseconds, 
%                           assumed to be the same for each image if t is a
%                           scalar
%       'colormap':         colormap
%       'spatialScaling':   spatial scaling factor 
%       'timeAcceleration': temporal scaling factor
%       'textSize':         size of font
%       'showTime':         print times over images
%       'showGrid':         overlay grid on images
%       'overlayText':      structure containing strings (s.str) to print
%                           at offset (s.offset.x, s.offset.y) over images


%% Parse Inputs

% initiate
p = inputParser;

% set defaults
default.t               = 500;   % milliseconds
default.imlimits        = [min(ims(:)) max(ims(:))];
default.filedir         = pwd;
default.filename        = 'ims';
default.colormap        = gray;
default.spatialScaling  = 1; 
default.timeAcceleration = 1;
default.textSize        = 12;
default.showTime        = false;
default.showProgressBar = false;
default.showGrid        = false;
default.overlayText     = struct([]);

% add arguments
addRequired(p,'ims',@isnumeric);
addParameter(p,'t',default.t,@isnumeric);
addParameter(p,'imlimits',default.imlimits,@isnumeric);
addParameter(p,'filedir',default.filedir,@ischar);
addParameter(p,'filename',default.filename,@ischar);
addParameter(p,'spatialScaling',default.spatialScaling,@isnumeric);
addParameter(p,'colormap',default.colormap,@isnumeric);
addParameter(p,'timeAcceleration',default.timeAcceleration,@isnumeric);
addParameter(p,'textSize',default.textSize,@isnumeric);
addParameter(p,'showTime',default.showTime,@islogical);
addParameter(p,'showProgressBar',default.showProgressBar,@islogical);
addParameter(p,'showGrid',default.showGrid,@islogical);
addParameter(p,'overlayText',default.overlayText,@isstruct);

% parse
parse(p,ims,varargin{:});
t = p.Results.t;
imlimits = p.Results.imlimits;
filedir = p.Results.filedir;
filename = p.Results.filename;
cmap = p.Results.colormap;
spatialScaling = p.Results.spatialScaling;
timeAcceleration = p.Results.timeAcceleration;
textSize = p.Results.textSize;
showTime = p.Results.showTime;
showProgressBar = p.Results.showProgressBar;
showGrid = p.Results.showGrid;
overlayText = p.Results.overlayText;

% further intialisation
if (length(t)==1)
    t = t * ((1:size(ims,3))-1); % in milliseconds
end


%% Setup


filepath = fullfile(filedir,[filename,'.gif']);

maxX = size(ims,2);
maxY = size(ims,1);

textOffset = 10;

textColor = 0.6*[1 1 1];

screenSize = get(0,'screensize');
screenSizeX = screenSize(3);
screenSizeY = screenSize(4)-100;
maxScaling = min(screenSizeX/maxX,screenSizeY/maxY);

magnification = min(spatialScaling,maxScaling)*100;

dt = diff(t);
dt(length(dt)+1) = dt(end);


%% Make Gif

hFig = figure('Visible','off');

hImg = imshow(ims(:,:,1),imlimits,'border','tight','initialMagnification',magnification);

if showProgressBar,
    hProgBar = add_progress_bar(ims,'showHistory');
end

if (~isempty(overlayText)),
    hTxt = text(round(overlayText(1).offset.x),...
                round(overlayText(1).offset.y),...
                overlayText(1).str,...
                'Color',textColor,...
                'HorizontalAlignment','left',...
                'FontSize',textSize);
end

colormap( cmap ),

nImage = size(ims,3);

for iImage = 1:nImage, 

    hImg.CData = ims(:,:,iImage);
        
    if (showGrid),
       overlay_gridlines( hImg.CData ); 
    end
    
    if exist('hProgBar','var'),
        update_progress_bar(hProgBar,iImage,nImage);
    end
    
    if (showTime)  % TODO: change time text to handle that gets string updated
        figText = sprintf('%.2f s',t(iImage)*1e-3);
        text(maxX-textOffset,maxY-textOffset,figText,'Color',textColor,'HorizontalAlignment','Right','FontSize',textSize)
    end
    
    if (~isempty(overlayText)),
         hTxt.String = overlayText(iImage).str;
    end
    
    pause(0.05) % short pause required so that subsequent call to getframe() 
                % references the current frame and not the previous

    frame = getframe(gca);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);

    delayTime = dt(iImage)/1000/timeAcceleration;
    
    if iImage ==1,
        imwrite(A,map,filepath,'gif','LoopCount',Inf,'DelayTime',delayTime);
    else
        imwrite(A,map,filepath,'gif','WriteMode','append','DelayTime',delayTime);
    end
    
end

close(hFig),

end  % ims2gif(...)


function overlay_gridlines( im )

M = size(im,1);
N = size(im,2);

step = 10;
xGrid = round(linspace(0,M,round(M/step))) + 0.5;
step = mode(diff(xGrid));
yGrid = round(linspace(0,N,round(N/step))) + 0.5;

hold on

for k = xGrid,
    x = [1 N];
    y = [k k];
    plot(x,y,'Color',0.2*[1,1,1],'LineStyle','-');
    plot(x,y,'Color','k','LineStyle',':');
end

for k = yGrid,
    x = [k k];
    y = [1 M];
    plot(x,y,'Color',0.2*[1,1,1],'LineStyle','-');
    plot(x,y,'Color','k','LineStyle',':');
end

hold off

end  % overlay_gridlines(...)


function hProgBar = add_progress_bar(ims,showHistory)


if ~exist('showHistory','var')
    showHistory = 'noHistory';
end

[maxY,maxX,nImage] = size(ims);

iImage = 1;

hProgBar.Outline  = rectangle('Position',[3,maxY-5,maxX-5,3],...
                              'Curvature',1,...
                              'LineWidth',1.5,...
                              'EdgeColor',[0.9 0.9 0.9],...
                              'FaceColor',[0.2 0.2 0.2]);
                          
if strcmp(showHistory,'showHistory')
    hProgBar.History  = rectangle('Position',hProgBar.Outline.Position,...
                              'Curvature',hProgBar.Outline.Curvature,...
                              'LineWidth',hProgBar.Outline.LineWidth,...
                              'EdgeColor',hProgBar.Outline.EdgeColor,...
                              'FaceColor',[0.4 0.4 0.4]);  
    hProgBar.History.Position(3) = hProgBar.History.Position(1);
end

hProgBar.Progress = rectangle('Position',hProgBar.Outline.Position,...
                              'Curvature',hProgBar.Outline.Curvature,...
                              'LineWidth',hProgBar.Outline.LineWidth,...
                              'EdgeColor',hProgBar.Outline.EdgeColor,...
                              'FaceColor',[0.7 0.7 0.7]); 
hProgBar.Progress.Position(3) = max( hProgBar.Outline.Position(4),...
                                     1/nImage * hProgBar.Outline.Position(3) );

update_progress_bar(hProgBar,iImage,nImage);
                      
end  % add_progress_bar(...)


function update_progress_bar(hProgBar,iImage,nImage)

widthOutlBar = hProgBar.Outline.Position(3);

if isfield(hProgBar,'History')    
    widthProgBar = max( hProgBar.Outline.Position(4),...
                        1/nImage * hProgBar.Outline.Position(3));
    hProgBar.History.Position(3) = ...
                            (iImage-1)/(nImage-1) * (widthOutlBar-widthProgBar) ...
                          + widthProgBar;                    
end

if isfield(hProgBar,'Progress')   
    widthProgBar = hProgBar.Progress.Position(3);
    hProgBar.Progress.Position(1) = ...
                            (iImage-1)/(nImage-1) * (widthOutlBar-widthProgBar) ...
                          + hProgBar.Outline.Position(1);                      
end

end  % update_progress_bar(...)