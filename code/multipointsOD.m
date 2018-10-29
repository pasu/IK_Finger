clear
clc

posInfo = {'dX','dY'};
S = load('smooth.mat',posInfo{:});
dX = getfield(S,'dX');
dY = getfield(S,'dY');

targetX = dX(1)
targetY = dY(1);
targetM = [targetX; targetY];
dPart = 0.4;

% Jacobian Inverse matrix
syms dRadian1 dRadian2 dRadian3 dDistance1 dDistance2 dDistance3
f1 =  cos(dRadian1 + dRadian2 + dRadian3) * dDistance3 ...
    + cos(dRadian1 + dRadian2) * dDistance2 ...
    + cos(dRadian1) * dDistance1;

f2 =  sin(dRadian1 + dRadian2 + dRadian3) * dDistance3 ...
    + sin(dRadian1 + dRadian2) * dDistance2 ...
    + sin(dRadian1) * dDistance1;

dScope1Max = pi/3;
dScope1Min = -pi/3;
dMean1 = (dScope1Max+dScope1Min)/2;

dScope2Max = 0;
dScope2Min = -2*pi/3;
dMean2 = (dScope2Max+dScope2Min)/2;

dScope3Max =0;
dScope3Min = -2*pi/3;
dMean3 = (dScope3Max+dScope3Min)/2;

%  constrainted equation
f3 = sin(dRadian2)*sin(dRadian2) + sin(dRadian3)*sin(dRadian3);

% f3 = ((dRadian1-dMean1)/(dMean1-dScope1Max))^2/3 ...
%             + ((dRadian2-dMean2)/(dMean2-dScope2Max))^2/3 ...
%             + ((dRadian3-dMean3)/(dMean3-dScope3Max))^2/3 ;
        
% f3 = 1/6 * ((dScope1Max-dScope1Min)/((dScope1Max-r1)*(r1-dScope1Min))  ...
%     + (dScope2Max-dScope2Min)/((dScope2Max-r2)*(r2-dScope2Min)) ...
%     + (dScope3Max-dScope3Min)/((dScope3Max-r3)*(r3-dScope3Min)));


f3J = jacobian(f3, [dRadian1,dRadian2,dRadian3]); 
% J
fxy = jacobian([f1,f2], [dRadian1,dRadian2,dRadian3]); 
fxyT = transpose(fxy);
fM = fxy * fxyT;
fxyV = inv(fM);
%fResult1 = fxyT * fxyV;
% J#:transpose(J) * inverse(J *  transpose(J))
fResult = fxyT / fM;
% I - J * J#
I = eye(3);
fResultOD = I - fResult * fxy;

% parameters
dRadian1 =0.1;
dRadian2 =0.1;
dRadian3 = 0.1;

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

% guess the start position
% dRadian1 = dGuessRadian1;
% dRadian2 = dGuessRadian2;
% dRadian3 = dGuessRadian3;
dRadian1 = 1;
 dRadian2 =  2*pi/3;
  dRadian3 = 2*pi/3;
clf;

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

rM = 0;
rP = -pi/3;
rD = -pi/3;
qR = [rM;rP; rD];
q0 = [double(dRadian1);double(dRadian2);double(dRadian3)];

for nPoint = 1:1:L
    targetX = dX(L-nPoint+1);% L-nPoint+1
    targetY = dY(L-nPoint+1);
    targetM = [targetX; targetY];
    
    clf;
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

for nTime = 1:1:50
    targetT = targetM - [Xq; Yq];

    fValue = vpa(subs(fResult));
    fValueOD = vpa(subs(fResultOD));

    qS = [dRadian1;dRadian2; dRadian3];
    part1 =  fValue * targetT;
    gValue = vpa(subs(f3J));
    part2 = fValueOD * transpose(gValue)  * (dPart);
    qS = qS + part1 + part2;

    dRadian1 = double((qS(1, 1)));
    dRadian2 = double((qS(2, 1)));
    dRadian3 = double((qS(3, 1)));
    
    if(dRadian1>dScope1Max)
        dRadian1 = rM;
    end
    if(dRadian1<dScope1Min)
        dRadian1 = rM;
    end

     if(dRadian2>dScope2Max)
        dRadian2 = rP;
    end
    if(dRadian2<dScope2Min)
        dRadian2 = rP;
    end
    
     if(dRadian3>dScope3Max)
        dRadian3 = rP;%0.5*dRadian2;
    end
    if(dRadian3<dScope3Min)
       dRadian3 = rP;%0.5*dRadian2;
    end
    
 %   dRadian1 = 0;
%    dRadian3 = double(-pi/2 - dRadian1-dRadian2);
%    dRadian3 = 2*dRadian2/3;
    

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

   
    distance = sqrt((targetX - Xq)*(targetX - Xq) + (targetY - Yq)*(targetY - Yq));
    if distance<1
        disp(distance);
        disp(nTime);
        break;
    end;
end
  disp('Press a key !')  % Press a key here.You can see the message 'Paused: Press any key' in        % the lower left corner of MATLAB window.
    pause;
end
%fValue = subs(fxy,[dRadian1,dRadian2,dRadian3,dDistance1,dDistance2,dDistance3],...
%    [pi/3.0,pi/4.0,pi/6.0,39.8,22.4,15.8]);

%fxyT = transpose(fValue);

%fxyV = inv(fM)