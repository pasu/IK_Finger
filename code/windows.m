function windows()
f=figure;
set(f,'WindowButtonDownFcn',@mytestcallback)
dDistance1 = 39.8;
dDistance2 = 22.4;
dDistance3 = 15.8;
set(gca,'XLim',[-(80)  dDistance1 + dDistance2 + dDistance3+10]);
set(gca,'YLim',[-( dDistance1 + dDistance2 + dDistance3) 65]);
grid on
     
function mytestcallback(hObject,~)
pos=get(gca,'CurrentPoint');
disp(['You clicked X:',num2str(pos(1)),', Y:',num2str(pos(2))]);