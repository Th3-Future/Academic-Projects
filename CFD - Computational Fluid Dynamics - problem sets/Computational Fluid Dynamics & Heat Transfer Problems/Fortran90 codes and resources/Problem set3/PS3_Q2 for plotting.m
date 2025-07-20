clear; clf;
soln=dlmread('PS3_Q2 t4.txt');
T=soln(1:80,2:81);
y=soln(1:80,1);
x=soln(81,2:81);

% clf; set(gcf,'position',[100 100 800 500]); 
mainDataAxes = axes('position',[0.1  0.1  0.8    0.85]);

% plot map
axes(mainDataAxes);
h=imagesc('XData',x(1:80),'YData',y(1:80),'CData',T);
% each data poit is exactly what it is. without interpolation.

% change colormap; to make the background white.
 j = jet; % colormap is simply an array of colors coded in RGB triplet.
 colormap(j); 
 
 g=colorbar;
 ylabel(g,'Temperature','FontSize',16,'Rotation',90)
 set( g, 'YDir', 'reverse' );
 % set the range of data that will be colorred
axHdl = get(h,'Parent'); % get the axes parent handle of the surface
get(axHdl,'CLim') % check out the range of colored data

% both x and y range from 1 to 68;
xlabel('x')
ylabel('y')
title(['Heat Map for  TOL_s_s= 5E-10  (@ time = 2561.89sec)'],'FontSize', 13)