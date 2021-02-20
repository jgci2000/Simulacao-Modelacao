%ellipse
% Function which draws an ellipse.
%
% LAST UPDATED by Andy French Jan 2013
%
% Syntax:  [x,y,r,theta,ecc,xf_r,xf_l,e_area] = ellipse(xc,yc,a,b,N)
%
% xc,yc    Cartesian coordinate of ellipse centre
% a        Distance from xc to ellipse furthest extent in x direction
% b        Distance from yc to ellipse furthest extent in y direction
% N        Number of points used to compute ellipse (equally spaced in
%          parameter phi used to parametrtize the ellipse cartesian
%          equation)
%
% x        1xN vector of Cartesian x coordinates
% y        1xN vector of Cartesian y coordinates
% r        1xN vector of polar r coordinates from focus at (xc + ecc*a,yc)
% theta    1xN vector of polar angles measured anticlockwise from a line parrallel
%          to the x axis x coordinates, passing through the focus
% ecc      Ellipse eccentricity ecc = sqrt( 1 - b^2 / a^2 )
% xf_r     x coordinate of right focus xc + ecc*a.
% xf_l     x coordinate of left focus xc - ecc*a.
% e_area   Area of ellipse = pi*a*b

function [x,y,r,theta,ecc,xf_r,xf_l,e_area] = ellipse(xc,yc,a,b,N)

%Define parametric angle in radians
phi = linspace(0,2*pi,N);

%Cinstruct x and y coordinates
x = xc + a*cos(phi);
y = yc + b*sin(phi);

%Eccentricity
ecc = sqrt( 1 - (b/a)^2 );

%Define polar coordinates, centre at right focus
r = sqrt( (x-xc-a*ecc).^2 + (y-yc).^2 );
theta = atan2( y-yc, x-xc - a*ecc );

%Test polar equation of ellipse
r_test = a*(1-ecc^2)./(1+ecc*cos(theta));
d = mean(abs(r-r_test));

%Foci and ellipse area
xf_r = xc + a*ecc;
xf_l = xc - a*ecc;
e_area = pi*a*b;


%%

%Copy and paste the contents of this function into the command window to
%test the ellipse function
function ellipse_test

xc = 3;
yc = -2;
a = 5;
b = 3;
N = 100;

[x,y,r,theta,ecc,xf_r,xf_l,e_area] = ellipse(xc,yc,a,b,N);

figure('name','ellipse.m test','color',[1 1 1])
plot(x,y,'b');
axis equal
view(2);
grid on
xlabel('x')
ylabel('y')
title({['Ellipse: a = ',num2str(a),', b = ',num2str(b),', xc = ',num2str(xc),', yc = ',num2str(yc)],...
    ['Eccentricity = ',num2str(ecc),', focus #1 = (',num2str(xf_l),',',num2str(yc),'), focus #2 = (',...
    num2str(xf_r),',',num2str(yc),'), area = ',num2str(e_area)]})
hold on
plot(xf_l,yc,'b*');
plot(xf_r,yc,'b*');
plot(xc,yc,'r*')
xlims = get(gca,'xlim');
ylims = get(gca,'ylim');
plot([0,0],ylims,'r-');
plot(xlims,[0,0],'r-');

%End of code
