function []=PlotPipeArrow(x,y,r,vfr,lw,rgb,tag)

alfa=atan2(y(2)-y(1),x(2)-x(1));

xp = mean(x);
yp = mean(y);

xa = [mean(x)-r*cos(alfa),mean(x)+r*cos(alfa)];
ya = [mean(y)-r*sin(alfa),mean(y)+r*sin(alfa)];

%color
R = rgb(1);
G = rgb(2);
B = rgb(3);

%start and end
plot([x(1),x(2)],[y(1),y(2)],'linewidth',lw,'color',[R,G,B],'tag',tag);
%plotting > or <
if(vfr>0)
    plot([xa(1)+r*cos(alfa-pi/2),xp],[ya(1)+r*sin(alfa-pi/2),yp],'linewidth',lw,'color',[R,G,B],'tag',tag); % \
    plot([xa(1)+r*cos(alfa+pi/2),xp],[ya(1)+r*sin(alfa+pi/2),yp],'linewidth',lw,'color',[R,G,B],'tag',tag); % /
elseif(vfr<0)
    plot([xa(2)+r*cos(alfa-pi/2),xp],[ya(2)+r*sin(alfa-pi/2),yp],'linewidth',lw,'color',[R,G,B],'tag',tag); % /
    plot([xa(2)+r*cos(alfa+pi/2),xp],[ya(2)+r*sin(alfa+pi/2),yp],'linewidth',lw,'color',[R,G,B],'tag',tag); % \
end

end
