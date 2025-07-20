clear; clf;
solnvr=dlmread('PS5_Q1vr.txt');
solnvz=dlmread('PS5_Q1vz.txt');
solnw=dlmread('PS5_Q1w.txt');
w=solnw(1:201,2:22);
Vr=solnvr(1:201,2:22);
Vz=solnvz(1:201,2:22);
z=solnvr(1:201,1);
r=solnvr(202,1:21);

figure(1)
axes('position',[0.1  0.1  0.4   0.85]);
sv = 12;                                                            % ‘Step’ Value
q=quiver(r, z(1:sv:end), Vr(1:sv:end,1:end), Vz(1:sv:end,1:end));
%'AutoScale','on', 'AutoScaleFactor', 2,'MaxHeadSize',10, 'LineStyle','-'
set(q,'AutoScale','on', 'AutoScaleFactor', 1, 'LineWidth',1,'MaxHeadSize',0.015)
axis tight
xlabel('r',"FontSize",13);
ylabel('z',"FontSize",13);
title(['Velocity Flow Field Vector'],'FontSize', 13);

figure(2)
% clf; set(gcf,'position',[100 100 800 500]); 
mainDataAxes = axes('position',[0.1  0.1  0.5   0.85]);

% plot map
axes(mainDataAxes);
axis tight
h=imagesc('XData',r(1:21),'YData',z(1:201),'CData',w);
% each data poit is exactly what it is. without interpolation.

% change colormap; to make the background white.
 j = jet; % colormap is simply an array of colors coded in RGB triplet.
 colormap(j); 
 
 g=colorbar;
 ylabel(g,'Vorticity','FontSize',16,'Rotation',90)
 set( g, 'YDir', 'reverse' );
 % set the range of data that will be colorred
axHdl = get(h,'Parent'); % get the axes parent handle of the surface
get(axHdl,'CLim') % check out the range of colored data

% both x and y range from 1 to 68;
xlabel('r',"FontSize",13);
ylabel('z',"FontSize",13);
title(['Vorticity Contour Field'],'FontSize', 13);