clear;

caseName = 'ferto';
% caseName = 'ujhermes';
% caseName = 'ky1';
% caseName = 'C-town';
% caseName = 'nagycenk';
% caseName = 'tomalom';
% caseName = 'lovo';
% caseName = 'balf';
% caseName = 'becsidomb';
% caseName = 'villasor';
% caseName = 'Anytown';

% plot settings
pumpRadius = 0.001;
presRadius = 0.001;
tankSize = 0.002;
valveSize = 0.00005;
lineWidth = 1.2;
nodeMarkerSize = 5;
backgroundColor = [0.9,0.9,0.9];
closedColor = [1.0,0.0,0.0];

% colorbar settings

%blackBody, blackBodyExt, cividis, coolWarmBent, coolWarmSmooth, inferno, jet, kindlmann, kindlmannExt, magma, plasma, viridis
%discrete: lines, prism
colorMapName = 'plasma'; 
colorBarText = 'Head [m]';
margin = 0.03; % margins around plot
colorDelta = 0.15; %space for colorbar
colorPos = "east"; %position of colorbar (east or south)
colorBarFontSize = 12;
colorElement = "Node"; % "Node" or "Pipe" or "All"

caseFolder = '../../Networks/Sopron/';
addpath('../../Plot');
projectFolder = pwd;

operatingSystem = ispc;
if(operatingSystem)
    slashSign = '\';
else
    slashSign = '/';
end

% loading the INP file
LoadINP;

% loading data for coloring
LoadColors;

fig=figure;
hold on;
% adjusting the size of the plot window
if(colorPos == "east")
    xlim([-margin,max(coordX)+margin+colorDelta]);
    ylim([-margin,max(coordY)+margin]);
elseif(colorPos == "south")
    xlim([-margin,max(coordX)+margin]);
    ylim([-margin-colorDelta,max(coordY)+margin]);
else
    disp('colorPos can be: south | east'); 
end
xPos = 50;
yPos = 50;
xyRatio = diff(xlim)/diff(ylim);
if(xyRatio<16/9)
   ySize = 1000;
   xSize = ySize*xyRatio;
else
   xSize = 1800;
   ySize = xSize/xyRatio;
end
set(gcf, 'Position',  [xPos, yPos, xSize, ySize]);
axis off; 
set(gcf,'InvertHardCopy','off','Color',backgroundColor);
%colorbar
cb = colorbar;
colorMap = importdata(['../../Plot/ColorMaps/',colorMapName,'.txt']);
%adjusting to discrete colormaps
if(colorMapName == "lines" || colorMapName == "prism")
    colorMap = repmat(colorMap,ceil(colorMapNumber/size(colorMap,1)),1);
    colorMap = colorMap(1:colorMapNumber,:);
    colorMap(:,1) = linspace(0,1,size(colorMap,1));
end
    % else
    colormap(colorMap(:,2:4));
    cb.FontSize = colorBarFontSize;
    cb.Label.String=colorBarText;
    for i=1:length(colorBarTicks)
       cb.TickLabels{i} = round(colorBarTicks(i),-round(log(max(abs(colorBarTicks)))/log(10))+2); 
    end
    if(colorPos == "south")
        cb.Orientation = 'horizontal';
        set(cb,'position',[0.1,0.08,0.8,0.025]);
    elseif(colorPos == "east")
        set(cb,'position',[0.9,0.1,0.025,0.8]);
    else
        disp('!!! WARNING !!! available colorPos: south | east'); 
    end
% end

if(colorElement == "Node")
   pipeData = zeros(size(pipe));
elseif(colorElement == "Pipe")
   nodeData = zeros(size(node));
end

% start of plotting
for i=1:pipeCounter
    if(colorElement == "Pipe" || colorElement == "All")
        if(~isnan(pipeColor(i)))
            r = min(1,max(0,interp1(colorMap(:,1),colorMap(:,2),pipeColor(i),'spline')));
            g = min(1,max(0,interp1(colorMap(:,1),colorMap(:,3),pipeColor(i),'spline')));
            b = min(1,max(0,interp1(colorMap(:,1),colorMap(:,4),pipeColor(i),'spline')));
        else
            r = closedColor(1);
            g = closedColor(2);
            b = closedColor(3);
        end
    else
        r = 0.2;
        g = 0.2;
        b = 0.2;
    end
    x = [node(pipe(i).nodeFromIdx).coordX,node(pipe(i).nodeToIdx).coordX];
    y = [node(pipe(i).nodeFromIdx).coordY,node(pipe(i).nodeToIdx).coordY];
    if(isnan(pipeData(i)))
        text = "OFF";
    else
        text = pipeData(i);
    end
    plot([x(1),(x(1)+x(2))/2,x(2)],[y(1),(y(1)+y(2))/2,y(2)],'k','linewidth',lineWidth,'color',[r,g,b],'tag',join([strrep(pipe(i).ID,'_','\_'),'  ',text]));
end
for i=1:pumpCounter
    if(colorElement == "All")
        if(~isnan(pumpColor(i)))
            r = min(1,max(0,interp1(colorMap(:,1),colorMap(:,2),pumpColor(i),'spline')));
            g = min(1,max(0,interp1(colorMap(:,1),colorMap(:,3),pumpColor(i),'spline')));
            b = min(1,max(0,interp1(colorMap(:,1),colorMap(:,4),pumpColor(i),'spline')));
        else
            r = closedColor(1);
            g = closedColor(2);
            b = closedColor(3);
        end
    else
        r = 0.5;
        g = 0.5;
        b = 0.5;
    end
    if(isnan(pumpData(i)))
        text = "OFF";
    else
        text = pumpData(i);
    end
    PlotPump([node(pump(i).nodeFromIdx).coordX,node(pump(i).nodeToIdx).coordX],[node(pump(i).nodeFromIdx).coordY,node(pump(i).nodeToIdx).coordY],pumpRadius,1,[r,g,b],join([strrep(pump(i).ID,'_','\_'),'  ',text]));
end
for i=1:presCounter
    if(colorElement == "All")
        if(~isnan(presColor(i)))
            r = min(1,max(0,interp1(colorMap(:,1),colorMap(:,2),presColor(i),'spline')));
            g = min(1,max(0,interp1(colorMap(:,1),colorMap(:,3),presColor(i),'spline')));
            b = min(1,max(0,interp1(colorMap(:,1),colorMap(:,4),presColor(i),'spline')));
        else
            r = closedColor(1);
            g = closedColor(2);
            b = closedColor(3);
        end
    else
        r = 0.5;
        g = 0.5;
        b = 0.5;
    end
    if(isnan(presData(i)))
        text = "OFF";
    else
        text = presData(i);
    end
    PlotPres(pres(i).coordX,pres(i).coordY,presRadius,1,90,[r,g,b],join([strrep(pres(i).ID,'_','\_'),'  ',text]));
end
for i=1:poolCounter
    if(colorElement == "All")
        if(~isnan(poolColor(i)))
            r = min(1,max(0,interp1(colorMap(:,1),colorMap(:,2),poolColor(i),'spline')));
            g = min(1,max(0,interp1(colorMap(:,1),colorMap(:,3),poolColor(i),'spline')));
            b = min(1,max(0,interp1(colorMap(:,1),colorMap(:,4),poolColor(i),'spline')));
        else
            r = closedColor(1);
            g = closedColor(2);
            b = closedColor(3); 
        end
    else
        r = 0.5;
        g = 0.5;
        b = 0.5;
    end
    if(isnan(poolData(i)))
        text = "OFF";
    else
        text = poolData(i);
    end
    PlotPool(pool(i).coordX,pool(i).coordY,tankSize,1,0,[r,g,b],join([strrep(pool(i).ID,'_','\_'),'  ',text]));
end
for i=1:valveCounter
    if(colorElement == "All")
        if(~isnan(valveColor(i)))
            r = min(1,max(0,interp1(colorMap(:,1),colorMap(:,2),valveColor(i),'spline')));
            g = min(1,max(0,interp1(colorMap(:,1),colorMap(:,3),valveColor(i),'spline')));
            b = min(1,max(0,interp1(colorMap(:,1),colorMap(:,4),valveColor(i),'spline')));
        else
            r = closedColor(1);
            g = closedColor(2);
            b = closedColor(3); 
        end
        if(isnan(valveData(i)))
            text = "OFF";
        else
            text = valveData(i);
        end
    else
        r = 0.5;
        g = 0.5;
        b = 0.5;
        text = "";
    end
    PlotValve([node(valve(i).nodeFromIdx).coordX,node(valve(i).nodeToIdx).coordX],[node(valve(i).nodeFromIdx).coordY,node(valve(i).nodeToIdx).coordY],valveSize,1,[r,g,b],join([strrep(valve(i).ID,'_','\_'),'  ',text]));
end
for i=1:nodeCounter
    if(colorElement == "Node" || colorElement == "All")
        if(~isnan(nodeColor(i)))
            r = min(1,max(0,interp1(colorMap(:,1),colorMap(:,2),nodeColor(i),'spline')));
            g = min(1,max(0,interp1(colorMap(:,1),colorMap(:,3),nodeColor(i),'spline')));
            b = min(1,max(0,interp1(colorMap(:,1),colorMap(:,4),nodeColor(i),'spline')));
        else
            r = closedColor(1);
            g = closedColor(2);
            b = closedColor(3);
        end
    else
        r = 0.2;
        g = 0.2;
        b = 0.2;
    end
    if(isnan(nodeData(i)))
        text = "OFF";
    else
        text = nodeData(i);
    end
    plot(node(i).coordX,node(i).coordY,'ko','linewidth',1,'markersize',nodeMarkerSize,'markerfacecolor',[r,g,b],'color',[0,0,0],'tag',join([strrep(node(i).ID,'_','\_'),'  ',text]));
end

datacursormode on
dcm = datacursormode(gcf);
set(dcm,'UpdateFcn',@ShowTag)

SaveTightFigure(fig,caseName);