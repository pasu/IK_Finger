clear
clc

syms dRadian1 dRadian2 dRadian3 dDistance1 dDistance2 dDistance3
f1 =  cos(dRadian1 + dRadian2 + dRadian3) * dDistance3 ...
    + cos(dRadian1 + dRadian2) * dDistance2 ...
    + cos(dRadian1) * dDistance1;

f2 =  sin(dRadian1 + dRadian2 + dRadian3) * dDistance3 ...
    + sin(dRadian1 + dRadian2) * dDistance2 ...
    + sin(dRadian1) * dDistance1;

dDistance1 = 39.8;
dDistance2 = 22.4;
dDistance3 = 15.8;

%position table for guess
dGapM = double(pi/30);
dGapP = double(pi/30);
dGapD = double(pi/20);

dDistanceSquare = 100000;
dGuessRadian1 = 0.1;
dGuessRadian2 = 0.1;
dGuessRadian3 = 0.1;

dIterationM = int32(2*pi/3.0 / dGapM);
dIterationP = int32(pi/3.0 / dGapP);
dIterationD = 1;

index = 1;
dMArray = zeros(1, dIterationM*dIterationP*dIterationD);
dPArray = zeros(1, dIterationM*dIterationP*dIterationD);
dDArray = zeros(1, dIterationM*dIterationP*dIterationD);
dX = zeros(1, dIterationM*dIterationP*dIterationD);
dY = zeros(1, dIterationM*dIterationP*dIterationD);

d1 = 39.8;
d2 = 22.4;
d3 = 15.8;

dRadian1 = 0.0;
dRadian2 = 0.01;
dRadian3 = 0.01;
color = [1.0,0.788,0.055];
for dRadian1 =  -pi/3.0:dGapM:pi/3.0;
    for dRadian2 =  -2.0*pi/3.0:dGapP:0.0;
        %for dRadian3 =  -2.0*pi/3.0:dGapD:0.0;
        dRadian3 = -pi/2 - dRadian1 - dRadian2;
        if dRadian3 > 0
            continue;
        end
        if dRadian3 < -2.0*pi/3.0
            continue;
        end
            dX(index)= vpa(subs(f1)) ;
            dY(index) = vpa(subs(f2));
            plot(dX(index), dY(index), 'ko', 'MarkerSize', 1);hold on;
            daspect([1 1 1]);
            set(gca,'XLim',[-(20) d1+d2+d3+10]);
            set(gca,'YLim',[-(d1+d2+d3) 65]);
            grid on
            
            index = index+1;
       % end
    end
end

mX = mean(dX);
mY = mean(dY);
colorspec = {[1.0,0.0,0.0]; [0.0,1.0,0.0]; [0.0,0.0,1.0]; [0.5,1.0,1.0]; [0.0,0.0,0.0]};
dGapM = pi/3.0;
dGapP = pi/4.0;
dGapD = pi/5.0;

%for radianM =  -pi/3.0:dGapM:pi/3.0;
radianM = pi/3.0;
    x0 = d1 * cos(radianM);
    y0 = d1 * sin(radianM);
    plot([0,x0],[0,y0],'r');hold on
    for radianP =  -2.0*pi/3.0:dGapP:0.0;
        x1 = x0 + d2 * cos(radianP+radianM);
        y1 = y0 + d2 * sin(radianP+radianM);
        plot([x0,x1],[y0,y1],'g');
        plot(x1,y1, 'go', 'MarkerSize', 10);hold on
        %for radianD =  -2.0*pi/3.0:dGapD:0.0;
            radianD =  -pi/2 - radianM - radianP;
            if radianD > 0
                continue;
            end
            if radianD < -2.0*pi/3.0
                continue;
            end
            x2 = x1 + d3 * cos(radianD + radianP+radianM);
            y2 = y1 + d3 * sin(radianD + radianP+radianM);
            plot(x2,y2, 'bo', 'MarkerSize', 10);
            plot([x1,x2],[y1,y2],'b');hold on
            daspect([1 1 1]);
            set(gca,'XLim',[-(20) d1+d2+d3+10]);
            set(gca,'YLim',[-(d1+d2+d3) 65]);
            xlabel('X (mm)');
            ylabel('Y (mm)');
            title('Points visited by the tip on rotation of  \theta_M, \theta_P and \theta_D,  \theta_M +  \theta_P +  \theta_D = -90');
            grid on
        %end
        %pause(0.001);
    end
%end

radianM = -pi/3.0;
    x0 = d1 * cos(radianM);
    y0 = d1 * sin(radianM);
    plot([0,x0],[0,y0],'r');hold on
    for radianP =  0:dGapP:0.0;
        x1 = x0 + d2 * cos(radianP+radianM);
        y1 = y0 + d2 * sin(radianP+radianM);
        plot([x0,x1],[y0,y1],'g');
        plot(x1,y1, 'go', 'MarkerSize', 10);hold on
        %for radianD =  -2.0*pi/3.0:dGapD:0.0;
        radianD =  -pi/2 - radianM - radianP;
        if radianD > 0
            continue;
        end
        if radianD < -2.0*pi/3.0
            continue;
        end
            x2 = x1 + d3 * cos(radianD + radianP+radianM);
            y2 = y1 + d3 * sin(radianD + radianP+radianM);
            plot(x2,y2, 'bo', 'MarkerSize', 10);
            plot([x1,x2],[y1,y2],'b');hold on
            daspect([1 1 1]);
            set(gca,'XLim',[-(20) d1+d2+d3+10]);
            set(gca,'YLim',[-(d1+d2+d3) 65]);
            xlabel('X (mm)');
            ylabel('Y (mm)');
            title('Points visited by the tip on rotation of  \theta_M, \theta_P and \theta_D,  \theta_M +  \theta_P +  \theta_D = -90');
            grid on
        %end
        %pause(0.001);
    end

plot(mX,mY,'-p','MarkerFaceColor','red', 'MarkerSize', 30);
    
