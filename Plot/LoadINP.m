
fileID = fopen([caseFolder, caseName,'.inp'],'r');

line = fgetl(fileID);
nodeCounter = 0;
presCounter = 0;
poolCounter = 0;
pipeCounter = 0;
pumpCounter = 0;
valveCounter = 0;
coordCounter = 0;

while(ischar(line))
   if(line == "[JUNCTIONS]")
      line = fgetl(fileID);
      while(~isempty(line) && line(1) == ";")%skipping over header
          line = fgetl(fileID);
      end
      while(line ~= "")
          nodeCounter = nodeCounter + 1;
          data = strsplit(line);
          data = data(data~=""); % clearing empty cells
          node(nodeCounter).ID = string(data(1));
          node(nodeCounter).elev = str2double(data(2));
          if(length(data)>3)
              node(nodeCounter).demand = str2double(data(3));
          else
              node(nodeCounter).demand = "";
          end
          if(length(data)>4)
              node(nodeCounter).pattern = string(data(4));
          else
              node(nodeCounter).pattern = "";
          end
          line = fgetl(fileID);
      end
   end
   if(line == "[RESERVOIRS]")
      line = fgetl(fileID);
      while(~isempty(line) && line(1) == ";")%skipping over header
          line = fgetl(fileID);
      end
      while(line ~= "")
          presCounter = presCounter + 1;
          data = strsplit(line);
          data = data(data~=""); % clearing empty cells
          pres(presCounter).ID = string(data(1));
          pres(presCounter).head = str2double(data(2));
          if(length(data)>3)
              pres(presCounter).pattern = string(data(3));
          else
              pres(presCounter).pattern = "";
          end
          pres(presCounter).type = "PressurePoint";
          % we are creating a node for the pressure point with the same ID
          nodeCounter = nodeCounter + 1;
          node(nodeCounter).ID = string(data(1));
          node(nodeCounter).elev = str2double(data(2));
          node(nodeCounter).demand = 0;
          node(nodeCounter).pattern = "";
          line = fgetl(fileID);
      end
   end
   if(line == "[TANKS]")
      line = fgetl(fileID);
      while(~isempty(line) && line(1) == ";")%skipping over header
          line = fgetl(fileID);
      end
      while(line ~= "")
          poolCounter = poolCounter + 1;
          data = strsplit(line);
          data = data(data~=""); % clearing empty cells
          pool(poolCounter).ID = string(data(1));
          pool(poolCounter).elev = str2double(data(2));
          pool(poolCounter).initLevel = str2double(data(3));
          pool(poolCounter).minLevel = str2double(data(4));
          pool(poolCounter).maxLevel = str2double(data(5));
          pool(poolCounter).diameter = str2double(data(6));
          if(length(data)>7)
              pool(poolCounter).minVol = string(data(7));
          else
              pool(poolCounter).minVol = "";
          end
          if(length(data)>8)
              pool(poolCounter).pattern = string(data(8));
          else
              pool(poolCounter).pattern = "";
          end
          pool(poolCounter).type = "Pool";
          % we are creating a node for the pool with the same ID
          nodeCounter = nodeCounter + 1;
          node(nodeCounter).ID = string(data(1));
          node(nodeCounter).elev = str2double(data(2));
          node(nodeCounter).demand = 0;
          node(nodeCounter).pattern = "";
          line = fgetl(fileID);
      end
   end
   if(line == "[PIPES]")
      line = fgetl(fileID);
      while(~isempty(line) && line(1) == ";")%skipping over header
          line = fgetl(fileID);
      end
      while(line ~= "")
          pipeCounter = pipeCounter + 1;
          data = strsplit(line);
          data = data(data~=""); % clearing empty cells
          pipe(pipeCounter).ID = string(data(1));
          pipe(pipeCounter).nodeFrom = string(data(2));
          pipe(pipeCounter).nodeTo = string(data(3));
          pipe(pipeCounter).length = str2double(data(4));
          pipe(pipeCounter).diameter = str2double(data(5));
          pipe(pipeCounter).roughness = str2double(data(6));
          if(length(data)>7)
              pipe(pipeCounter).minorLoss = str2double(data(7));
          else
              pipe(pipeCounter).minorLoss = "";
          end
          if(length(data)>8)
              pipe(pipeCounter).status = str2double(data(8));
          else
              pipe(pipeCounter).minorLoss = "";
          end
          pipe(pipeCounter).type = "Pipe";
          line = fgetl(fileID);
      end
   end
   if(line == "[PUMPS]")
      line = fgetl(fileID);
      while(~isempty(line) && line(1) == ";")%skipping over header
          line = fgetl(fileID);
      end
      while(line ~= "")
          pumpCounter = pumpCounter + 1;
          data = strsplit(line);
          data = data(data~=""); % clearing empty cells
          pump(pumpCounter).ID = string(data(1));
          pump(pumpCounter).nodeFrom = string(data(2));
          pump(pumpCounter).nodeTo = string(data(3));
          pump(pumpCounter).pumpType = string(data(4));
          pump(pumpCounter).parameter = str2double(data(5));
          pump(pumpCounter).type = "Pump";
          line = fgetl(fileID);
      end
   end
   if(line == "[VALVES]")
      line = fgetl(fileID);
      while(~isempty(line) && line(1) == ";")%skipping over header
          line = fgetl(fileID);
      end
      while(line ~= "")
          valveCounter = valveCounter + 1;
          data = strsplit(line);
          data = data(data~=""); % clearing empty cells
          valve(valveCounter).ID = string(data(1));
          valve(valveCounter).nodeFrom = string(data(2));
          valve(valveCounter).nodeTo = string(data(3));
          valve(valveCounter).diameter = str2double(data(4));
          valve(valveCounter).valveType = string(data(5));
          valve(valveCounter).setting = str2double(data(6));
          valve(valveCounter).minorLoss = str2double(data(7));
          valve(valveCounter).type = "Valve";
          line = fgetl(fileID);
      end
   end
   if(line == "[COORDINATES]")
      line = fgetl(fileID);
      while(~isempty(line) && line(1) == ";")%skipping over header
          line = fgetl(fileID);
      end
      while(line ~= "")
          coordCounter = coordCounter + 1;
          data = strsplit(line);
          data = data(data~=""); % clearing empty cells
          coordName(coordCounter) = string(data(1));
          if(caseFolder == "../../Networks/Sopron/")
              coordX(coordCounter) = str2double(data(3));
              coordY(coordCounter) = -str2double(data(2));
          else              
              coordX(coordCounter) = str2double(data(2));
              coordY(coordCounter) = -str2double(data(3));
          end
          line = fgetl(fileID);
      end
   end
   line = fgetl(fileID);
end

% finding node indicies for pipes
for i=1:length(node) % this is for searching node IDs efficiently
    nodeID(i) = node(i).ID;
end
for i=1:presCounter
   pres(i).nodeIdx = find(nodeID==pres(i).ID);
end
for i=1:poolCounter
   pool(i).nodeIdx = find(nodeID==pool(i).ID);
end
for i=1:pipeCounter
   pipe(i).nodeFromIdx = find(nodeID==pipe(i).nodeFrom);
   pipe(i).nodeToIdx = find(nodeID==pipe(i).nodeTo);
end
for i=1:pumpCounter
   pump(i).nodeFromIdx = find(nodeID==pump(i).nodeFrom);
   pump(i).nodeToIdx = find(nodeID==pump(i).nodeTo);
end
for i=1:valveCounter
   valve(i).nodeFromIdx = find(nodeID==valve(i).nodeFrom);
   valve(i).nodeToIdx = find(nodeID==valve(i).nodeTo);
end


% normalizing x-y coordinates to [0,1] without distortion
delta = max([max(coordX)-min(coordX),max(coordY)-min(coordY)]);
coordX = abs(coordX - min(coordX))/delta;
coordY = abs(max(coordY) - coordY)/delta; % Flipping the plot

% finding node x-y coordinates
for i=1:nodeCounter
    idx = find(node(i).ID==coordName);
    node(i).coordX = coordX(idx);
    node(i).coordY = coordY(idx);
end
for i=1:presCounter
    idx = find(pres(i).ID==coordName);
    pres(i).coordX = coordX(idx);
    pres(i).coordY = coordY(idx);
end
for i=1:poolCounter
    idx = find(pool(i).ID==coordName);
    pool(i).coordX = coordX(idx);
    pool(i).coordY = coordY(idx);
end

fclose(fileID);


