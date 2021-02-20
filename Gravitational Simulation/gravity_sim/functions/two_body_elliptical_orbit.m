%two_body_elliptical_orbit
% Function which generates all the parameters to calculate mutually elliptical orbits
% caused by Newtonian gravitation given masses m1,m2 and maximum separation
% s.
%
% NOTE: Initial x,y,z coordinates are in AU and velocities are in AU per
% Year. Masses are in Solar masses.
%
% Syntax: [ux,vx1,vx2,...
%          uy,vy1,vy2,...
%          uz,vz1,vz2,...
%          xc,x1,x2,...
%          yc,y1,y2,...
%          zc,z1,z2,...
%          theta,theta_dot,t,j,J,E,P] =...
%    two_body_elliptical_orbit(  G, m1, m2, s, ecc, theta0, h0, hdot, d, alpha, N )
%
% For Newtonian Gravity G = 39.43 AU^3 M-1 Yr^-2  where M is the mass of
% the Sun. If another set of length (L), mass (M) and time (T) units are used, G must
% be modified accordingly with units of L^3 M^-1 T^-2. For a galactic
% simulation T = 3.2e15s, L = 4.8e20 m, M = 6e42 kg (i.e. parameters of the
% Milky Way). In this case G = 37.
%
% s corresponds to the maximum separation of masses
%
% ecc is the eccentricity of the orbit. 0 < ecc < 1.
% ecc = 0 corresponds to a circular orbit.
%
% theta0 is the initial phase angle of the orbit /radians, measured anti-clockwise
% from the semi-major axis.
%
% hdot is the velocity of the centre of mass of the system. This is a constant of the motion
% as is defined as an [ux,uy,uz] vector.
%
% h0 is the initial coordinate of the centre of mass of the system [h0x,h0y,h0z]
%
% d = [dx,dy,dz] is the direction of the semi-major axis of the ellipse.
%
% alpha is the anti-clockwise rotation of the ellipse about the direction vector /radians
%
% N is the number of data points used to compute the orbit. The first and
% last coordinate will be be repeated to close the loop.
%
% Outputs:
%
% [ux,uy,uz]        Velocities of centre of mass
% [vx1,vy1,vz1]     Velocities of mass #1
% [vx2,vy2,vz2]     Velocities of mass #2
% [xc,yc,zc]        Position vector of centre of mass
% [x1,y1,z1]        Position vector of mass #1
% [x2,y2,z2]        Position vector of mass #2
% theta             Elliptical angle, starting with theta_0 /radians
% theta_dot         Rate of change of elliptical angle /radians /time unit T
% t                 Time of orbit /time unit T
% j                 Reduced angular momentum of masses (ignoring rotation
%                   of centre of mass about coordinate system origin).
% J                 Total angular momentum of system
% E                 Total energy of the system
% P                 Orbital period /time unit T
%
% LAST UPDATED by A. French April 2013

function [ux,vx1,vx2,...
    uy,vy1,vy2,...
    uz,vz1,vz2,...
    xc,x1,x2,...
    yc,y1,y2,...
    zc,z1,z2,...
    theta,theta_dot,t,j,J,E,P] =...
    two_body_elliptical_orbit(  G, m1, m2, s, ecc, theta0, h0, hdot, d, alpha,N )

%Semi-major axis a of ellipse traced out by mass separation vector
a = s/(1+ecc);

%Ellipse polar angle
theta = linspace( theta0, theta0 + 2*pi, N );

%Time /Yr
t = sqrt( a^3 * (1 - ecc^2)^3 / ( G *(m1+m2) )) * theta_integral( theta, theta0, ecc );

%Position of centre of masses
xc = h0(1) + t.*hdot(1) ;
yc = h0(2) + t.*hdot(2) ;
zc = h0(3) + t.*hdot(3) ;

%Velocity of centre of masses
ux = hdot(1);
uy = hdot(2);
uz = hdot(3);

%Mass separation vector components
w = a*(1-ecc^2)./(1 + ecc*cos(theta) );
wx = w.*cos(theta);
wy = w.*sin(theta);
wz = zeros(1,N);

%Separation vector rate of change
wdot = sqrt( G*(m1+m2)/(a*(1-ecc^2) )) * (1 +ecc*cos(theta) );
wdotx = wdot.*( cos(theta)*ecc.*sin(theta)./(1+ecc*cos(theta))  - sin(theta) );
wdoty = wdot.*( sin(theta).*ecc.*sin(theta)./(1+ecc*cos(theta))  + cos(theta) );
wdotz = zeros(1,N);

%Rotate to desired orientation
[wx,wy,wz] = point(wx,wy,wz,1,0,0,d(1),d(2),d(3));
[wx,wy,wz] = vrot(wx,wy,wz,d(1),d(2),d(3),d(1),d(2),d(3),alpha);
[wdotx,wdoty,wdotz] = point(wdotx,wdoty,wdotz,1,0,0,d(1),d(2),d(3));
[wdotx,wdoty,wdotz] = vrot(wdotx,wdoty,wdotz,d(1),d(2),d(3),d(1),d(2),d(3),alpha);

%Position of masses
x1 = xc - m2*wx/(m1+m2);
y1 = yc - m2*wy/(m1+m2);
z1 = zc - m2*wz/(m1+m2);
x2 = xc + m1*wx/(m1+m2);
y2 = yc + m1*wy/(m1+m2);
z2 = zc + m1*wz/(m1+m2);

%Velocity of masses
vx1 = hdot(1) - m2*wdotx/(m1+m2);
vy1 = hdot(2) - m2*wdoty/(m1+m2);
vz1 = hdot(3) - m2*wdotz/(m1+m2);
vx2 = hdot(1) + m1*wdotx/(m1+m2);
vy2 = hdot(2) + m1*wdoty/(m1+m2);
vz2 = hdot(3) + m1*wdotz/(m1+m2);

%Rate of change of polar angle
theta_dot = ( (1 + ecc*cos(theta) ).^2 ) * sqrt( G*(m1+m2)/( a^3 * (1-ecc^2)^3 ) );

%Angular momentum
j = ( m1*m2/(m1+m2) ) * cross( [wx(1),wy(1),wz(1)], [wdotx(1),wdoty(1),wdotz(1)] );
J = j + (m1+m2) * cross( h0,hdot );

%Total energy
E = 0.5*( m2 + m2 ) * dot(hdot,hdot) - G*m1*m2/(2*a);

%Orbital period P
P = sqrt( 4*pi^2 * a^3 / (G*(m1+m2)) );

%%

%Numerical integration of the time, theta integral for the orbit
function tau = theta_integral( theta, theta0, ecc )
if length(theta)==1
    tau = 0;
else
    [Y_pp,tau] = calculus( theta,(1 + ecc*cos(theta)).^-2,theta,'integ','spline',0);
    [Y_pp,tau0] = calculus( theta,(1 + ecc*cos(theta)).^-2,theta0,'integ','spline',0);
    tau = tau-tau0;
end

%%

%CALCULUS
% Function that fits a piecwise polynomial to an array of [X,Y] data and
% thus computes the derivative and cumulative integral of Y(X) at all
% points X. The function is based on the code written on pages 292 and 294
% of "Mastering Matlab 6" by Duane Hanselman & Bruce Littlefield.
%
% LAST UPDATED by Andrew French. 21/09/2004.
%
% Syntax: [Y_pp,Y_calc]=calculus(X,Y,X0,type,poly,C)
%
% X        - Vector of data X variable
% Y        - Vector of data Y(X) variable
% X0       - Vector of X variables corresponding to that of the output.
% type     - 'integ' or 'diff'. Choses integration or diffrentiation of Y(X).
% poly     - 'spline' or 'hermite' . Choses type of piecewise polynomial
%            computed. Hermite polynomials tend to be better for
%            non-smooth data and tends to be 'shape preserving.'
% C        - If type=='integ', this paramater gives the initial value of the
%            integral corresponding to X(1).
%
% Y_pp     - Piecewise polynomial fitted to Y(X) evaluated at X0.
% Y_calc   - Vector of the same dimensions of X0, giving the
%            differential, or integral of Y(X) evaluated at all values of X0.
%
% Example: X=linspace(-2,2,20); Y=X.^2;
% [Y_pp,Y_int]=calculus(X,Y,X,'integ','spline',0);
% [Y_int,YY]=calculus(X,Y_int,X,'diff','spline',0);
% plot(X,Y); hold on;
% plot(X,YY,'g');
% plot(X,Y_pp,'ro');
% legend('Y=X^2','Integral then derivative of Y=X^2','PP fit of Y=X^2');
% xlabel('X'); ylabel('Y');

function [Y_pp,Y_calc]=calculus(X,Y,X0,type,poly,C)

%Compute piecewise polynomial. Cubic spline is the default
if strcmp(poly,'hermite')==1
    pp=pchip(X,Y);
else
    pp=spline(X,Y);
end

%Check if integration or differentiation is required. Differentiation is
%the default option.
if strcmp(type,'integ')==1
    
    %Check that C is a scalar
    if prod(size(C))~=1
        error('C must be a scalar!')
    end
    
    %Take apart piecewise polynomial
    [br,co,npy,nco]=unmkpp(pp);
    
    %Scale factors for integration
    sf=nco:-1:1;
    
    %Integral coefficients
    ico=[co./sf(ones(npy,1),:) zeros(npy,1)];
    
    %Integral spline has higher order
    nco=nco+1;
    
    %Integration constant
    ico(1,nco)=C;
    
    %Find constant terms in polynomials
    for k=2:npy
        ico(k,nco)=polyval(ico(k-1,:),br(k)-br(k-1));
    end
    
    %Build pp form for integral
    ppi=mkpp(br,ico);
    
    %Evaluate integral at values of X0
    Y_calc=ppval(ppi,X0);
    
    %Evaluate piecewise polynomial
    Y_pp=ppval(pp,X0);
else
    
    %Take apart piecewise polynomial
    [br,co,npy,nco]=unmkpp(pp);
    
    %Scale factors for differentiation
    sf=nco-1:-1:1;
    
    %Derivative coefficients
    dco=sf(ones(npy,1),:).*co(:,1:nco-1);
    
    %Build pp form for derivative
    ppd=mkpp(br,dco);
    
    %Evaluate derivative at values of X0
    Y_calc=ppval(ppd,X0);
    
    %Evaluate piecewise polynomial
    Y_pp=ppval(pp,X0);
end

%End of code