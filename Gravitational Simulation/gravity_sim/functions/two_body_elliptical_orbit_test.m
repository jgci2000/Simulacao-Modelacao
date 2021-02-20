function two_body_elliptical_orbit_test

%Define inputs
G = 39.43;
m1 = 1;
m2 = 2; 
s = 4;
ecc = 0.7;
theta0 = pi/4;
h0 = [0,0,0];
hdot = [0,0,1];
d = [0,1,0];
alpha = pi/30;
N = 500;

%Generate orbits
[ux,vx1,vx2,...
          uy,vy1,vy2,...
          uz,vz1,vz2,...
          xc,x1,x2,...
          yc,y1,y2,...
          zc,z1,z2,...
          theta,theta_dot,t,j,J,E,P] =...
    two_body_elliptical_orbit(  G, m1, m2, s, ecc, theta0, h0, hdot, d, alpha,N );

%Generate title string showing orbital parameters
title_str = {['m1 = ',num2str(m1),', m2 = ',num2str(m2),', s = ',num2str(s),', ecc = ',num2str(ecc),', period = ',num2str(P),', G = ',num2str(G)],...
    ['theta0 = ',num2str(theta0),', h0 = [',num2str(h0),'], hdot = [',num2str(hdot),'], d = [',num2str(d),'], alpha = ',num2str(alpha)],...
    ['J = [',num2str(J),'], E = ',num2str(E)]};

%Plot orbits
figure('color',[1 1 1],'name','two body elliptical orbit test',...
    'units','normalized')
plot3(x1,y1,z1,'r'); hold on;
plot3(x2,y2,z2,'b');
plot3(xc,yc,zc,'k');
plot3(x1(1),y1(1),z1(1),'ro');
plot3(x2(1),y2(1),z2(1),'bo');
plot3(xc(1),yc(1),zc(1),'ko');
plot3(x1(N),y1(N),z1(N),'r*');
plot3(x2(N),y2(N),z2(N),'b*');
plot3(xc(N),yc(N),zc(N),'k*');
sc = 0.1;
quiver3(x1(1),y1(1),z1(1),sc*vx1(1),sc*vy1(1),sc*vz1(1),0,'r');
quiver3(x1(N/2),y1(N/2),z1(N/2),sc*vx1(N/2),sc*vy1(N/2),sc*vz1(N/2),0,'r');
quiver3(x2(1),y2(1),z2(1),sc*vx2(1),sc*vy2(1),sc*vz2(1),0,'b');
quiver3(x2(N/2),y2(N/2),z2(N/2),sc*vx2(N/2),sc*vy2(N/2),sc*vz2(N/2),0,'b');
xlim([-s/2,s/2]);
ylim([-s/2,s/2]);
zlim([-s/2,s/2]);
axis vis3d;
grid on
xlabel('x')
ylabel('y')
zlabel('z')
title(title_str)
set(gca,'units','normalized','outerposition',[0,0,1,1]);
set(gca,'position',[0.0,0.13,1,0.7])
print( gcf, ['Elliptical orbit test plot',...
        '. ',strrep(datestr(now),':','-'),'.png'],'-dpng','-r300' )

%Plot theta vs time
figure('color',[1 1 1],'name','two body elliptical orbit test time',...
    'units','normalized')
plot(theta,t);
xlabel('Polar angle /radians')
ylabel('Orbit time /years')
title(title_str)
set(gca,'units','normalized','outerposition',[0,0,1,1]);
set(gca,'position',[0.12,0.15,0.8,0.7])
print( gcf, 'elliptical orbit test time','-dpng','-r300' )
close(gcf);

%End of code