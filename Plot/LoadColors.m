% checking the folder exist
if(~exist([projectFolder,slashSign,caseName]))
    disp(['Folder (',caseName,') not found']);
end

% loading the data for coloring
if(exist([projectFolder,slashSign,caseName,slashSign,'Node.txt']))
%     nodeColor = importdata([projectFolder,slashSign,caseName,slashSign,'Node.txt']);
    fileID = fopen([projectFolder,slashSign,caseName,slashSign,'Node.txt'],'r');
    nodeColor = fscanf(fileID,'%f');
    fclose(fileID);
    nodeData = nodeColor;
else
    nodeColor = 0;
    nodeData = 0;
end
if(exist([projectFolder,slashSign,caseName,slashSign,'Pipe.txt']))
%     pipeColor = importdata([projectFolder,slashSign,caseName,slashSign,'Pipe.txt']);
    fileID = fopen([projectFolder,slashSign,caseName,slashSign,'Pipe.txt'],'r');
    pipeColor = fscanf(fileID,'%f');
    fclose(fileID);
    pipeData = pipeColor;
else
    pipeColor = 0;
    pipeData = 0;
end
if(exist([projectFolder,slashSign,caseName,slashSign,'Pool.txt']))
%     poolColor = importdata([projectFolder,slashSign,caseName,slashSign,'Pool.txt']);
    fileID = fopen([projectFolder,slashSign,caseName,slashSign,'Pool.txt'],'r');
    poolColor = fscanf(fileID,'%f');
    fclose(fileID);
    poolData = poolColor;
else
    poolColor = 0;
    poolData = 0;
end
if(exist([projectFolder,slashSign,caseName,slashSign,'Pres.txt']))
%     presColor = importdata([projectFolder,slashSign,caseName,slashSign,'Pres.txt']);
    fileID = fopen([projectFolder,slashSign,caseName,slashSign,'Pres.txt'],'r');
    presColor = fscanf(fileID,'%f');
    fclose(fileID);
    presData = presColor;
else
    presColor = 0;
    presData = 0;
end
if(exist([projectFolder,slashSign,caseName,slashSign,'Pump.txt']))
%     pumpColor = importdata([projectFolder,slashSign,'',caseName,slashSign,'Pump.txt']);
    fileID = fopen([projectFolder,slashSign,caseName,slashSign,'Pump.txt'],'r');
    pumpColor = fscanf(fileID,'%f');
    fclose(fileID);
    pumpData = pumpColor;
else
    pumpColor = 0;
    pumpData = 0;
end

allColor = [nodeColor;pipeColor;presColor;poolColor;pumpColor];

if(colorMapName == "lines" || colorMapName == "prism")
   colorMapNumber = length(unique(allColor));
end

% defining the labels for colorbar
if(colorElement == "Pipe")
    if(max(pipeColor)-min(pipeColor) == 0)
        colorBarTicks = repmat(max(pipeColor),11,1);
    else
        colorBarTicks = min(pipeColor):(max(pipeColor)-min(pipeColor))/10:max(pipeColor);
    end
elseif(colorElement == "Node")
    if(max(nodeColor)-min(nodeColor) == 0)
        colorBarTicks = repmat(max(nodeColor),11,1);
    else
        colorBarTicks = min(nodeColor):(max(nodeColor)-min(nodeColor))/10:max(nodeColor);
    end
elseif(colorElement == "All")
    if(max(allColor)-min(allColor) == 0)
        colorBarTicks = repmat(max(allColor),11,1);
    else
        colorBarTicks = min(allColor):(max(allColor)-min(allColor))/10:max(allColor);
    end
elseif(colorElement == "Both")
    colorBarTicks = 0:0.1:1;
else
   disp("!!! WARNING !!! Available colorELement: Pipe, Node, All"); 
   colorBarTicks = [0,1];
end
% reset to [0,1] interval
if(colorElement == "All")
    if((max(allColor)-min(allColor)) == 0)
        denom = max(allColor);
    else
        denom = max(allColor)-min(allColor);
    end
    nodeColor = (nodeColor-min(allColor))/denom;
    pipeColor = (pipeColor-min(allColor))/denom;
    poolColor = (poolColor-min(allColor))/denom;
    presColor = (presColor-min(allColor))/denom;
    pumpColor = (pumpColor-min(allColor))/denom;
end
if(colorElement == "Node" || colorElement == "Both")
    if((max(nodeColor)-min(nodeColor)) == 0)
        denom = max(nodeColor);
    else
        denom = max(nodeColor)-min(nodeColor);
    end
    nodeColor = (nodeColor-min(nodeColor))/denom;
end
if(colorElement == "Pipe" || colorElement == "Both")
    if((max(pipeColor)-min(pipeColor)) == 0)
        denom = max(pipeColor);
    else
        denom = max(pipeColor)-min(pipeColor);
    end
    pipeColor = (pipeColor-min(pipeColor))/denom;
end





