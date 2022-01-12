% checking the folder exist
if(~exist([projectFolder,slashSign,'Network Data',slashSign,caseName]))
    disp(['Folder (','Network Data',slashSign,caseName,') not found']);
end

% loading the data for coloring
if(exist([projectFolder,slashSign,'Network Data',slashSign,caseName,slashSign,'Node.txt']))
    fileID = fopen([projectFolder,slashSign,'Network Data',slashSign,caseName,slashSign,'Node.txt'],'r');
    nodeColor = fscanf(fileID,'%f');
    fclose(fileID);
    nodeData = nodeColor;
else
    nodeColor = 0;
    nodeData = 0;
end
if(exist([projectFolder,slashSign,'Network Data',slashSign,caseName,slashSign,'Pipe.txt']))
    fileID = fopen([projectFolder,slashSign,'Network Data',slashSign,caseName,slashSign,'Pipe.txt'],'r');
    pipeColor = fscanf(fileID,'%f');
    fclose(fileID);
    pipeData = pipeColor;
else
    pipeColor = 0;
    pipeData = 0;
end
if(exist([projectFolder,slashSign,'Network Data',slashSign,caseName,slashSign,'Pool.txt']))
    fileID = fopen([projectFolder,slashSign,'Network Data',slashSign,caseName,slashSign,'Pool.txt'],'r');
    poolColor = fscanf(fileID,'%f');
    fclose(fileID);
    poolData = poolColor;
else
    poolColor = 0.2*ones(10,1);
    poolData = zeros(10,1);
end
if(exist([projectFolder,slashSign,'Network Data',slashSign,caseName,slashSign,'Pres.txt']))
    fileID = fopen([projectFolder,slashSign,'Network Data',slashSign,caseName,slashSign,'Pres.txt'],'r');
    presColor = fscanf(fileID,'%f');
    fclose(fileID);
    presData = presColor;
else
    presColor = 0.2*ones(10,1);
    presData = zeros(10,1);
end
if(exist([projectFolder,slashSign,'Network Data',slashSign,caseName,slashSign,'Pump.txt']))
    fileID = fopen([projectFolder,slashSign,'Network Data',slashSign,caseName,slashSign,'Pump.txt'],'r');
    pumpColor = fscanf(fileID,'%f');
    fclose(fileID);
    pumpData = pumpColor;
else
    pumpColor = 0.2*ones(10,1);
    pumpData = zeros(20,1);
end
if(exist([projectFolder,slashSign,'Network Data',slashSign,caseName,slashSign,'Valve.txt']))
    fileID = fopen([projectFolder,slashSign,'Network Data',slashSign,caseName,slashSign,'Valve.txt'],'r');
    valveColor = fscanf(fileID,'%f');
    fclose(fileID);
    valveData = valveColor;
else
    valveColor = 0.5*ones(5000,1);
    valveData = zeros(5000,1);
end

allColor = [nodeColor;pipeColor;presColor;poolColor;pumpColor;valveColor];

if(colorMapName == "lines" || colorMapName == "prism")
   colorMapNumber = length(unique(allColor));
end

% defining the labels for colorbar
if(~isLogColorMap)
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
    elseif(colorElement == "Valve")
        if(max(valveColor)-min(valveColor) == 0)
            colorBarTicks = repmat(max(valveColor),11,1);
        else
            colorBarTicks = min(valveColor):(max(valveColor)-min(valveColor))/10:max(valveColor);
        end
    elseif(colorElement == "All")
        if(max(allColor)-min(allColor) == 0)
            colorBarTicks = repmat(max(allColor),11,1);
        else
            colorBarTicks = min(allColor):(max(allColor)-min(allColor))/10:max(allColor);
        end
    elseif(colorElement == "Both")
        colorBarTicks = 0:0.1:1;
    elseif(colorElement == "None")
        colorBarTicks = 0:0.1:1;
    else
       disp("!!! WARNING !!! Available colorELement: Pipe, Node, All, None"); 
       colorBarTicks = [0,1];
    end
else
   colorBarTicks =  zeros(1,colorBarTicksNumber);
   lin = min(nodeColor):(max(nodeColor)-min(nodeColor))/(colorBarTicksNumber-1):max(nodeColor);
   colorBarTicks(1) = lin(1);
   colorBarTicks(end) = lin(end);
   for i=length(lin)-1:-1:2
      colorBarTicks(i) = (colorBarTicks(i+1)+colorBarTicks(1))/2; 
   end
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
if(colorElement == "Valve")
    if((max(valveColor)-min(valveColor)) == 0)
        denom = max(valveColor);
    else
        denom = max(valveColor)-min(valveColor);
    end
    valveColor = (valveColor-min(valveColor))/denom;
end



