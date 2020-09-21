clear;
% close all

% caseName = 'buk';
% caseName = 'vashegy';
% caseName = 'ky1';
% caseName = 'agyagosszergeny';
% caseName = 'csapod';
% caseName = 'lovo';

% caseName = 'pozsonyiut';

caseName = 'tomalom';

% caseName = 'villasor';
% caseName = 'kohegy';
% caseName = 'und';
% caseName = 'rojtokmuzsaj';

% caseName = 'Anytown';
% caseName = 'ky2';
% caseName = 'C-town';
% caseName = 'grid_16';

caseFolder = '../../Networks/Sopron/';

% plot settings
pumpRadius = 0.000015;
pumpColor = 0.8*ones(11,1);
pumpLineWidth = 1.6;
presRadius = 0.013;
presColor = 0.2*ones(10,1);
presAngle = [-180,0,90,90,0];
presLineWidth = 1.6;
poolSize = 0.026;
poolAngle = [0,-210,0,0,0,-75,0];
poolColor = 0.2*ones(10,1);
poolLineWidth = 1.6;
valveSize = 0.000022;
valveColor = 0.5*ones(2000,1);
valveLineWidth = 1.2;
lineWidth = 1.6;
nodeMarkerSize = 5;
backgroundColor = [1.0,1.0,1.0];
closedColor = [1.0,0.0,0.0];

% colorbar settings

%blackBody, blackBodyExt, cividis, coolWarmBent, coolWarmSmooth, inferno, jet, kindlmann, kindlmannExt, magma, plasma, viridis
%discrete: lines, prism
colorMapName = 'plasma'; 
colorBarText = 'Edge sensitvity [i]';
margin = 0.03; % margins around plot
colorDelta = 0.15; %space for colorbar
colorPos = "east"; %position of colorbar (east or south)
colorBarFontSize = 14;
colorElement = "Node"; % "Node" or "Pipe" or "All" or "None"
isLogColorMap = false; % WORKS ONLY WITH NODE COLORELEMENT
colorBarTicksNumber = 10;
colorBarVisible = 'on';

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
elseif(isLogColorMap)
    colorMap = [colorMap(1:colorBarTicksNumber-1:end,:);colorMap(end,:)];
    for i=2:size(colorMap,1)-1
       colorMap(i,1) = (colorMap(1,1) + colorMap(i+1,1))/2;
    end
end

colormap(colorMap(:,2:4));
cb.FontSize = colorBarFontSize;
cb.Label.String = colorBarText;
cb.TickLabels = round(colorBarTicks,-round(log(max(abs(colorBarTicks)))/log(10))+2); 
set(cb,'YTick',0:1/(colorBarTicksNumber-1):1)
cb.Visible = colorBarVisible;
if(colorPos == "south")
    cb.Orientation = 'horizontal';
    set(cb,'position',[0.1,0.08,0.8,0.025]);
elseif(colorPos == "east")
    set(cb,'position',[0.9,0.1,0.025,0.8]);
else
    disp('!!! WARNING !!! available colorPos: south | east'); 
end

% colorbar off
if(colorElement == "Node")
   pipeData = zeros(size(pipe));
elseif(colorElement == "Pipe")
   nodeData = zeros(size(node));
elseif(colorElement == "None")
    pipeData = zeros(size(pipe));
    nodeData = zeros(size(node));
%     pumpData = zeros(size(pump));
end

% start of plotting
x = zeros(length(pipe),3);
y = zeros(length(pipe),3);
for i=1:length(pipe)
    xx = [node(pipe(i).nodeFromIdx).coordX,node(pipe(i).nodeToIdx).coordX];
    x(i,:) = [xx(1),(xx(1)+xx(2))/2,xx(2)];
    yy = [node(pipe(i).nodeFromIdx).coordY,node(pipe(i).nodeToIdx).coordY];
    y(i,:) = [yy(1),(yy(1)+yy(2))/2,yy(2)];
end

plotObj = plot(x',y','k');

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
    elseif(colorElement == "None")
        r = 0.0;
        g = 0.0;
        b = 0.0;
    else
        r = 0.2;
        g = 0.2;
        b = 0.2;
    end
    if(isnan(pipeData(i)))
        text = "OFF";
    else
        text = pipeData(i);
    end
    plotObj(i).Color = [r,g,b];
    plotObj(i).Tag   = join([strrep(pipe(i).ID,'_','\_'),'  ',text]);
    plotObj(i).LineWidth = lineWidth;
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
    elseif(colorElement == "None")
        r = min(1,max(0,interp1(colorMap(:,1),colorMap(:,2),pumpColor(i),'spline')));
        g = min(1,max(0,interp1(colorMap(:,1),colorMap(:,3),pumpColor(i),'spline')));
        b = min(1,max(0,interp1(colorMap(:,1),colorMap(:,4),pumpColor(i),'spline')));
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
    PlotPump([node(pump(i).nodeFromIdx).coordX,node(pump(i).nodeToIdx).coordX],[node(pump(i).nodeFromIdx).coordY,node(pump(i).nodeToIdx).coordY],pumpRadius,pumpLineWidth,[r,g,b],join([strrep(pump(i).ID,'_','\_'),'  ',text]));
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
    elseif(colorElement == "None")
        r = min(1,max(0,interp1(colorMap(:,1),colorMap(:,2),valveColor(i),'spline')));
        g = min(1,max(0,interp1(colorMap(:,1),colorMap(:,3),valveColor(i),'spline')));
        b = min(1,max(0,interp1(colorMap(:,1),colorMap(:,4),valveColor(i),'spline')));
        text = "";
    else
        r = 0.5;
        g = 0.5;
        b = 0.5;
        text = "";
    end
    PlotValve([node(valve(i).nodeFromIdx).coordX,node(valve(i).nodeToIdx).coordX],[node(valve(i).nodeFromIdx).coordY,node(valve(i).nodeToIdx).coordY],valveSize,valveLineWidth,[r,g,b],join([strrep(valve(i).ID,'_','\_'),'  ',text]));
end

x = zeros(size(node));
y = zeros(size(node));
for i=1:length(node)
   x(i) = node(i).coordX; 
   y(i) = node(i).coordY; 
end
plotObj = plot([x;x],[y;y],'ko');

for i=1:nodeCounter
    if(colorElement == "Node" || colorElement == "All")
        if(~isnan(nodeColor(i)))
            r = min(1,max(0,interp1(colorMap(:,1),colorMap(:,2),nodeColor(i))));
            g = min(1,max(0,interp1(colorMap(:,1),colorMap(:,3),nodeColor(i))));
            b = min(1,max(0,interp1(colorMap(:,1),colorMap(:,4),nodeColor(i))));
        else
            r = closedColor(1);
            g = closedColor(2);
            b = closedColor(3);
        end
    elseif(colorElement == "None")
        r = 0.4;
        g = 0.4;
        b = 0.4;
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
    plotObj(i).MarkerFaceColor = [r,g,b];
    plotObj(i).Tag             = join([strrep(node(i).ID,'_','\_'),'  ',text]);
    plotObj(i).MarkerSize      = nodeMarkerSize;
    
%     plot(node(i).coordX,node(i).coordY,'o','linewidth',1,'markersize',nodeMarkerSize,'markerfacecolor',[r,g,b],'color',[0,0,0],'tag',join([strrep(node(i).ID,'_','\_'),'  ',text]));
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
        if(isnan(presData(i)))
            text = "OFF";
        else
            text = presData(i);
        end
    elseif(colorElement == "None")
        r = min(1,max(0,interp1(colorMap(:,1),colorMap(:,2),presColor(i),'spline')));
        g = min(1,max(0,interp1(colorMap(:,1),colorMap(:,3),presColor(i),'spline')));
        b = min(1,max(0,interp1(colorMap(:,1),colorMap(:,4),presColor(i),'spline')));
        text = "";
    else
        r = 0.5;
        g = 0.5;
        b = 0.5;
        text = "";
    end
    PlotPres(pres(i).coordX,pres(i).coordY,presRadius,presLineWidth,presAngle(i),[r,g,b],join([strrep(pres(i).ID,'_','\_'),'  ',text]));
end
for i=1:poolCounter
% for i=1:2
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
        if(isnan(poolData(i)))
            text = "OFF";
        else
            text = poolData(i);
        end
    elseif(colorElement == "None")
        r = min(1,max(0,interp1(colorMap(:,1),colorMap(:,2),poolColor(i),'spline')));
        g = min(1,max(0,interp1(colorMap(:,1),colorMap(:,3),poolColor(i),'spline')));
        b = min(1,max(0,interp1(colorMap(:,1),colorMap(:,4),poolColor(i),'spline')));
        text = "";
    else
        r = 0.5;
        g = 0.5;
        b = 0.5;
        text = "";
    end

    PlotPool(pool(i).coordX,pool(i).coordY,poolSize,poolLineWidth,poolAngle(i),[r,g,b],join([strrep(pool(i).ID,'_','\_'),'  ',text]));
end

datacursormode on
dcm = datacursormode(gcf);
set(dcm,'UpdateFcn',@ShowTag)

SaveTightFigure(fig,['Plots',slashSign,caseName]);