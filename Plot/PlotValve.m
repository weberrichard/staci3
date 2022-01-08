function []=PlotValve(x,y,r,lw,rgb,tag)

d=sqrt((x(2)-x(1))^2+(y(2)-y(1))^2);
alfa=atan2(y(2)-y(1),x(2)-x(1));
w=0.5*r; %width of valve
t=0.5*r; %size of T
tt = t/2;

xp = [mean(x)-r/2*cos(alfa),mean(x)+r/2*cos(alfa)];
yp = [mean(y)-r/2*sin(alfa),mean(y)+r/2*sin(alfa)];

%color
R = rgb(1);
G = rgb(2);
B = rgb(3);

%start and end
plot([x(1),xp(1)],[y(1),yp(1)],'linewidth',lw,'color',[R,G,B],'tag',tag);
plot([x(2),xp(2)],[y(2),yp(2)],'linewidth',lw,'color',[R,G,B],'tag',tag);
%plotting |><|
plot([xp(1)+w*cos(alfa-pi/2),xp(1)+w*cos(alfa+pi/2)],[yp(1)+w*sin(alfa-pi/2),yp(1)+w*sin(alfa+pi/2)],'linewidth',lw,'color',[R,G,B],'tag',tag); % -|
plot([xp(2)+w*cos(alfa-pi/2),xp(2)+w*cos(alfa+pi/2)],[yp(2)+w*sin(alfa-pi/2),yp(2)+w*sin(alfa+pi/2)],'linewidth',lw,'color',[R,G,B],'tag',tag); %    |-
plot([xp(1)+w*cos(alfa-pi/2),xp(2)+w*cos(alfa+pi/2)],[yp(1)+w*sin(alfa-pi/2),yp(2)+w*sin(alfa+pi/2)],'linewidth',lw,'color',[R,G,B],'tag',tag); %   \
plot([xp(2)+w*cos(alfa-pi/2),xp(1)+w*cos(alfa+pi/2)],[yp(2)+w*sin(alfa-pi/2),yp(1)+w*sin(alfa+pi/2)],'linewidth',lw,'color',[R,G,B],'tag',tag); %   /
%plotting T
plot([mean(x),mean(x)+t*cos(alfa+pi/2)],[mean(y),mean(y)+t*sin(alfa+pi/2)],'linewidth',lw,'color',[R,G,B],'tag',tag);
plot([mean(x)+t*cos(alfa+pi/2)+tt*cos(alfa+pi),mean(x)+t*cos(alfa+pi/2)+tt*cos(alfa)],[mean(y)+t*sin(alfa+pi/2)+tt*sin(alfa+pi),mean(y)+t*sin(alfa+pi/2)+tt*sin(alfa)],'linewidth',lw,'color',[R,G,B],'tag',tag);

end
