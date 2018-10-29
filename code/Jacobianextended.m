clear
clc

posInfo = {'dX','dY','dMArray','dPArray','dDArray'};
S = load('pos.mat',posInfo{:});
dX = getfield(S,'dX');
dY = getfield(S,'dY');
dMArray = getfield(S,'dMArray');
dPArray = getfield(S,'dPArray');
dDArray = getfield(S,'dDArray');

d1 = 39.8;
d2 = 22.4;
d3 = 15.8;

dPart = 0;
% the position
targetX = 25.4  ;
targetY = -20.7;
targetX = dX(123);
targetY = dY(123);
targetX = 5.4  ;
targetY = -20.7;
targetM = [targetX; targetY;0];

% Jacobian Inverse matrix
syms r1 r2 r3 d1 d2 d3
f1 =  cos(r1 + r2 + r3) * d3 ...
    + cos(r1 + r2) * d2 ...
    + cos(r1) * d1;

f2 =  sin(r1 + r2 + r3) * d3 ...
    + sin(r1 + r2) * d2 ...
    + sin(r1) * d1;

dScope1Max = pi/3;
dScope1Min = -pi/3;
dMean1 = (dScope1Max+dScope1Min)/2;

dScope2Max = 0;
dScope2Min = -2*pi/3;
dMean2 = (dScope2Max+dScope2Min)/2;

dScope3Max =0;
dScope3Min = -2*pi/3;
dMean3 = (dScope3Max+dScope3Min)/2;

%  three constrainted equation
% f3 = sin(r2)*sin(r2) + sin(r3)*sin(r3);
f3 = ((r1-dMean1)/(dMean1-dScope1Max))^2/3 ...
            + ((r2-dMean2)/(dMean2-dScope2Max))^2/3 ...
            + ((r3-dMean3)/(dMean3-dScope3Max))^2/3 ;

% f3 = 1/6 * ((dScope1Max-dScope1Min)/((dScope1Max-r1)*(r1-dScope1Min))  ...
%     + (dScope2Max-dScope2Min)/((dScope2Max-r2)*(r2-dScope2Min)) ...
%     + (dScope3Max-dScope3Min)/((dScope3Max-r3)*(r3-dScope3Min)));

fxy = jacobian([f1,f2], [r1,r2,r3]);

det1 =  det([fxy(3),fxy(5);fxy(4),fxy(6)]);
det2 = - det([fxy(1),fxy(5);fxy(2),fxy(6)]);
det3 =  det([fxy(1),fxy(3);fxy(2),fxy(4)]);

f3D = jacobian(f3, [r1,r2,r3]); 
G =  f3D(1) * det1 +  f3D(2) * det2 + f3D(3) * det3;
f3D = jacobian(G, [r1, r2,r3]);
% fxy(1) = -d1*sin(r1)-d2*sin(r1+r2)-d3*sin(r1+r2+r3);
% fxy(2) = -d2*sin(r1+r2)-d3*sin(r1+r2+r3);
% fxy(3) = -d3*sin(r1+r2+r3);
% fxy(4) = d1*cos(r1)+d2*cos(r1+r2)+d3*cos(r1+r2+r3);
% fxy(5) = d2*cos(r1+r2)+d3*cos(r1+r2+r3);
% fxy(6) = d3*cos(r1+r2+r3);
% f3D(1) = d1*d3*(sin(r2)^2-cos(r2)^2)*sin(r2+r3) - d1*d3*sin(r2)*cos(r2)*cos(r2+r3)+d2*d3*(sin(r2)^2-cos(r2)^2)*sin(r3)+d1*d2*cos(r2)*sin(r3)*cos(r3)+d1*d2*sin(r3)*cos(r3)*cos(r2+r3);
% f3D(2) = -d1*d2*sin(r2)*cos(r2)*cos(r2+r3) - d2*d3*sin(r2)*cos(r2)*cos(r3) -d1*d2*sin(r2)*(sin(r3)^2-cos(r3)^2) + d1*d2*sin(r3)*cos(r3)*cos(r2+r3) - d1*d2*(sin(r3)^2-cos(r3)^2)*sin(r2+r3);
% J
fxyE = ... 
[ fxy(1),fxy(3),fxy(5);...
fxy(2),fxy(4),fxy(6);...
f3D(1), f3D(2),f3D(3)];

% fxyE = ... 
% [ fxy(1),fxy(2),fxy(3);...
% fxy(4),fxy(5),fxy(6);...
% 0, f3D(1),f3D(2)];
fxyV = inv((fxyE));
% parameters
r1 =0;
r2 = 0;
r3 = 0.1;

d1 = 39.8;
d2 = 22.4;
d3 = 15.8;



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
% r1 = dGuessRadian1;
% r2 = dGuessRadian2;
% r3 = dGuessRadian3;

Xq = vpa(subs(f1));
Yq = vpa(subs(f2));

clf;
 x0 = d1 * cos(r1);
 y0 = d1 * sin(r1);
plot([0,x0],[0,y0],'r');hold on

x1 = x0 + d2 * cos(r1+r2);
y1 = y0 + d2 * sin(r1+r2);
plot([x0,x1],[y0,y1],'g');hold on

 x2 = x1 + d3 * cos(r1+r2+r3);
y2 = y1 + d3 * sin(r1+r2+r3);
plot([x1,x2],[y1,y2],'b');
%text(double(Xq),double(Yq),'start pos');
 hold on
 
% iteration research start:
plot(targetX, targetY, '*', 'MarkerSize', 8);

rM = 0;
rP = -pi/3;
rD = -pi/3;
qR = [rM;rP; rD];
q0 = [double(r1);double(r2);double(r3)];

for nTime = 1:1:20
    targetT = targetM - [Xq; Yq;0];
    
     fValue = vpa(subs(fxyV));

%     qS = [r1;r2; r3];
%     part1 =  fValue * targetT * dPart;
%     qS = qS + part1;

    qS = [r1;r2; r3];
    qS = qS +  fValue * targetT;

    r1 = double((qS(1, 1)));
    r2 = double((qS(2, 1)));
    r3 = double((qS(3, 1)));
    
    if(r1>dScope1Max)
        r1 = dScope1Max-0.1;
    end
    if(r1<dScope1Min)
        r1 = dScope1Min+0.1;
    end

     if(r2>dScope2Max)
        r2 = dScope2Max-0.1;
    end
    if(r2<dScope2Min)
        r2 = dScope2Min+0.1;
    end
%      r3 = double(-pi/2 - r1-r2);
     if(r3>dScope3Max)
        r3 = dScope3Max-0.1;%0.5*dRadian2;
    end
    if(r3<dScope3Min)
       r3 = dScope3Min+0.1;%0.5*dRadian2;
    end
%      r1 = 0;
%     r3 = double(-pi/2 - r1-r2);
%     r3 = 2*r2/3;

    Xq = vpa(subs(f1));
    Yq = vpa(subs(f2));
    
   x0 = d1 * cos(r1);
    y0 = d1 * sin(r1);
    plot([0,x0],[0,y0],'r');hold on

    x1 = x0 + d2 * cos(r1+r2);
    y1 = y0 + d2 * sin(r1+r2);
    plot([x0,x1],[y0,y1],'g');hold on

    x2 = x1 + d3 * cos(r1+r2+r3);
    y2 = y1 + d3 * sin(r1+r2+r3);
    plot([x1,x2],[y1,y2],'b');    
    text(double(Xq),double(Yq),num2str(nTime));
     set(gca,'XLim',[-(80)  d1 + d2 + d3+10]);
     set(gca,'YLim',[-( d1 + d2 + d3) 65]);
     grid on
     
     distance = double(sqrt((targetX - Xq)*(targetX - Xq) + (targetY - Yq)*(targetY - Yq)));
     resultD = [nTime r1/3.1415926*180 r2/3.1415926*180 r3/3.1415926*180 targetX targetY distance];
        s = num2str(resultD);
        disp(s);
        
    disp('Press a key !')  % Press a key here.You can see the message 'Paused: Press any key' in        % the lower left corner of MATLAB window.
    pause;

   
    if distance<1
         resultD = [r1/3.1415926*180 r2/3.1415926*180 r3/3.1415926*180 targetX targetY nTime];
        s = num2str(resultD);
        disp(s);
        break;
    end;
end