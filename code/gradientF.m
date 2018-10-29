clear
clc

dScope1Max = ((3.1415926/3));
dScope1Min =(( -3.1415926/3));
dMean1 = (( (dScope1Max+dScope1Min)/2));

% y1 = power(((dRadian1-dMean1)/(dMean1-dScope1Max)),2);
% y2 = ((dScope1Max-dScope1Min)/((dScope1Max-dRadian1).*(dRadian1-dScope1Min)))/2;

x = [-1:0.01:1];
y1 =  power(((x-dMean1)/(dMean1-dScope1Max)),2);
plot(x, y1);hold on;
y2 =  ((dScope1Max-dScope1Min)./((dScope1Max-x).*(x-dScope1Min)))./2;
plot(x, y2)
xlabel('\theta');
ylabel('H(\theta)');
title('Joint limit avoidance potential function');

% for x = -1:0.1:1
%     dRadian1 = x;
%     yV1 = vpa(subs(y1));
%     yV2 = vpa(subs(y2));
%     plot(dRadian1,y1,'go',dRadian1,y2,'b--o');
%     %plot(dRadian1,yV2);
% end
