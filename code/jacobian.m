clear
clc

posInfo = {'dX','dY','dMArray','dPArray','dDArray'};
S = load('pos.mat',posInfo{:});
dX = getfield(S,'dX');
dY = getfield(S,'dY');
dMArray = getfield(S,'dMArray');
dPArray = getfield(S,'dPArray');
dDArray = getfield(S,'dDArray');
% the position
targetX = 25;%dDistance1 + dDistance2 + dDistance3;
targetY = 46;
targetM = [targetX; targetY];

% Jacobian Inverse matrix
syms dRadian1 dRadian2 dRadian3 dDistance1 dDistance2 dDistance3
f1 =  cos(dRadian1 + dRadian2 + dRadian3) * dDistance3 ...
    + cos(dRadian1 + dRadian2) * dDistance2 ...
    + cos(dRadian1) * dDistance1;

f2 =  sin(dRadian1 + dRadian2 + dRadian3) * dDistance3 ...
    + sin(dRadian1 + dRadian2) * dDistance2 ...
    + sin(dRadian1) * dDistance1;

fxy = jacobian([f1,f2], [dRadian1,dRadian2,dRadian3]);
fxyT = transpose(fxy);
fM = fxy * fxyT;
fxyV = inv(fM);
%fResult1 = fxyT * fxyV;
fResult = fxyT / fM;
% parameters
dRadian1 = 0.01;
dRadian2 = 0.01;
dRadian3 = 0.01;

dDistance1 = 39.8;
dDistance2 = 22.4;
dDistance3 = 15.8;



%position table for guess
dGapM = double(pi/9);
dGapP = double(pi/9);
dGapD = double(pi/9);

dDistanceSquare = 100000;
dGuessRadian1 = 0.1;
dGuessRadian2 = 0.1;
dGuessRadian3 = 0.1;

L = length(dX);

for index = 1:1:L;
            Xq = dX(index);
            Yq = dY(index);
             dDistanceToTarget = (Xq-targetX) * (Xq-targetX) + (Yq-targetY) * (Yq-targetY);
                if dDistanceToTarget<dDistanceSquare
                    dGuessRadian1 =  dMArray(index);
                    dGuessRadian2 = dPArray(index);
                    dGuessRadian3 = dDArray(index);
                    dDistanceSquare = dDistanceToTarget;
                end
end

% guess the start position
dRadian1 = dGuessRadian1;
dRadian2 = dGuessRadian2;
dRadian3 = dGuessRadian3;

Xq = vpa(subs(f1));
Yq = vpa(subs(f2));

 x0 = dDistance1 * cos(dRadian1);
 y0 = dDistance1 * sin(dRadian1);
plot([0,x0],[0,y0],'r');hold on

x1 = x0 + dDistance2 * cos(dRadian1+dRadian2);
y1 = y0 + dDistance2 * sin(dRadian1+dRadian2);
plot([x0,x1],[y0,y1],'g');hold on

 x2 = x1 + dDistance3 * cos(dRadian1+dRadian2+dRadian3);
y2 = y1 + dDistance3 * sin(dRadian1+dRadian2+dRadian3);
plot([x1,x2],[y1,y2],'b');
%text(double(Xq),double(Yq),'start pos');
 hold on
 
% iteration research start:
plot(targetX, targetY, '*', 'MarkerSize', 8);

nIterationSize = 10;

dScope1Max = pi/3;
dScope1Min = -pi/3;

dScope2Max = 0;
dScope2Min = -2*pi/3;

dScope3Max =0;
dScope3Min = -2*pi/3;

for nTime = 1:1:10
    targetT = targetM - [Xq; Yq];

    fValue = vpa(subs(fResult));

    qS = [dRadian1;dRadian2; dRadian3];
    qS = qS + fValue * targetT;

    dRadian1 = double((qS(1, 1)));
    dRadian2 = double((qS(2, 1)));
    dRadian3 = double((qS(3, 1)));
    
%     if(dRadian1>dScope1Max)
%         dRadian1 = dScope1Max;
%     end
%     if(dRadian1<dScope1Min)
%         dRadian1 = dScope1Min;
%     end
%     
%      if(dRadian2>dScope2Max)
%         dRadian2 = dScope2Max;
%     end
%     if(dRadian2<dScope2Min)
%         dRadian2 = dScope2Min;
%     end
%     
%      if(dRadian3>dScope3Max)
%         dRadian3 = dScope3Min;
%     end
%     if(dRadian3<dScope3Min)
%         dRadian3 = dScope3Min;
%     end
    

    Xq = vpa(subs(f1));
    Yq = vpa(subs(f2));
    
   x0 = dDistance1 * cos(dRadian1);
    y0 = dDistance1 * sin(dRadian1);
    plot([0,x0],[0,y0],'r');hold on

    x1 = x0 + dDistance2 * cos(dRadian1+dRadian2);
    y1 = y0 + dDistance2 * sin(dRadian1+dRadian2);
    plot([x0,x1],[y0,y1],'g');hold on

    x2 = x1 + dDistance3 * cos(dRadian1+dRadian2+dRadian3);
    y2 = y1 + dDistance3 * sin(dRadian1+dRadian2+dRadian3);
    plot([x1,x2],[y1,y2],'b');    
    text(double(Xq),double(Yq),num2str(nTime));
     set(gca,'XLim',[-(80)  dDistance1 + dDistance2 + dDistance3+10]);
     set(gca,'YLim',[-( dDistance1 + dDistance2 + dDistance3) 65]);
     grid on
    disp('Press a key !')  % Press a key here.You can see the message 'Paused: Press any key' in        % the lower left corner of MATLAB window.
    pause;

   
    distance = sqrt((targetX - Xq)*(targetX - Xq) + (targetY - Yq)*(targetY - Yq))
    if distance<1
        break;
    end;
end
%fValue = subs(fxy,[dRadian1,dRadian2,dRadian3,dDistance1,dDistance2,dDistance3],...
%    [pi/3.0,pi/4.0,pi/6.0,39.8,22.4,15.8]);

%fxyT = transpose(fValue);

%fxyV = inv(fM)