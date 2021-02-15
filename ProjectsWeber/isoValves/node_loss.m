clear;
close all;

%blackBody, blackBodyExt, cividis, coolWarmBent, coolWarmSmooth, inferno, jet, kindlmann, kindlmannExt, magma, plasma, viridis
%discrete: lines, prism
% colorMapName = 'grayscale'; 
colorMapName = 'plasma';
colorMap = importdata(['../../Plot/ColorMaps/',colorMapName,'.col']);

 cases= ["villasor","ferto","sanchegy","buk","lovo","nagycenk","vashegy","varis","becsidomb","tomalom",...
    "szakov","kohegy","harka","pozsonyiut","sopronkovesd","dudlesz","ivan","agyagosszergeny","kofejto","simasag",...
    "acsad","csaford","nagylozs","balf","csapod","und","rojtokmuzsaj","brennberg","pusztacsalad","kutyahegy",...
    "nyarliget","meszlen","fertoujlak","gorbehalom","tozeggyarmajor","ebergoc","csillahegy","jerevan","gloriette",...
    "ohermes","ujhermes"];

% cases= ["ferto","sanchegy","lovo","nagycenk","vashegy","varis","becsidomb","tomalom",...
%     "szakov","kohegy","harka","pozsonyiut","sopronkovesd","dudlesz","ivan","agyagosszergeny","kofejto","simasag",...
%     "acsad","csaford","nagylozs","balf","csapod","und","rojtokmuzsaj","brennberg","pusztacsalad","kutyahegy",...
%     "nyarliget","meszlen","fertoujlak","gorbehalom","tozeggyarmajor","ebergoc","csillahegy","jerevan","gloriette",...
%     "ohermes","ujhermes"];

idx = 1:27;

cor_beta = zeros(size(idx));

k=1;
for i=idx
    beta = importdata(join(['Network Data/',cases(i),'/demand_loss_orig.txt'],''));
    gamma = importdata(join(['Network Data/',cases(i),'/vulner_orig.txt'],''));
    ev = importdata(join(['Network Data/',cases(i),'/segment_edge_orig.txt'],'')) + 1;
    max_node = max(max(ev));
    node_names = string(1:max_node);
    G = graph(ev(:,1),ev(:,2));
    G.Nodes.Name = node_names';
    pi = importdata(join(['Network Data/',cases(i),'/input_segment_orig.txt'],'')) + 1;
    pi = unique(pi);

   	Li = importdata(join(['Network Data/',cases(i),'/absolute_segment_length_orig.txt'],''));
    li = Li/sum(Li);
    Di = importdata(join(['Network Data/',cases(i),'/segment_demand_orig.txt'],''));
    di = Di/sum(Di);
    
    delta_1 = zeros(size(beta));
    delta_2 = zeros(size(beta));
    delta_3 = zeros(size(beta));
    % removing nodes one by one
    for j=1:max_node
        G2 = rmnode(G,j);
%         figure(1);
%         plot(G); title('G');
%         figure(2);
%         plot(G2); title('G2');
        c = conncomp(G2);
        [ispres,~] = find(string(pi')==G2.Nodes.Name);
        cc = unique(c(ispres));
        if(isempty(ispres))
            delta_1(j) = 1;
        else
            % simple node number
            delta_1(j) = (max_node-sum(sum(c'==cc)))/max_node;

            % relative pipe length
            v = li([1:j-1 j+1:end]);
            id = boolean(sum(c'==cc,2));
            vv = sum(sum(v(id)));
            delta_2(j) = 1-vv;
            
            % relative demand
            v = di([1:j-1 j+1:end]);
            id = boolean(sum(c'==cc,2));
            vv = sum(sum(v(id)));
            delta_3(j) = 1-vv;
        end
    end
    
    delta_1 = delta_1(beta<0.75);
    delta_2 = delta_2(beta<0.75);
    delta_3 = delta_3(beta<0.75);
    beta = beta(beta<0.75);

    cor_beta_1(k) = corr(delta_1,beta);
    cor_beta_2(k) = corr(delta_2,beta);
    cor_beta_3(k) = corr(delta_3,beta);
    
    k=k+1;
end

% figure();
% plot(delta,beta,'x');
% 
figure();
boxplot([cor_beta_1',cor_beta_2',cor_beta_3']);

