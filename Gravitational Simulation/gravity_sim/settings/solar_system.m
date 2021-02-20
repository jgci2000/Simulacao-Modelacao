%Solar System
%Two suns with equal mass. Each is surrounded by concentric rings which
%co-rotate with the circular orbit of the masses.

%Astronomical parameters in SI units
AU = 149.6e9;
G = 6.67e-11;
M_sun = 2e30;
Yr = 365*24*3600;

%Dimensionless number which controls the dynamics, and results from the
%scaling of mass, distance and time parameters to make them dimensionless.
G = G*Yr^2*M_sun/(AU^3);

%Strength of gravity
G_factor = 1;
G = G_factor*G;

%Scale the solar mass
M_Sun = 1;

%Set dimensions of space
Lmax = 40;

%Timestep
dt = 0.05;

%Gravitation exponent
Ng = 2;

%Repulsion exponent
Nr = 5;

%Repulsion magnitude
RM = 0.01;

%Boundary coefficient of restitution (if [] then there is
%no boundary)
k = [];

%Initial view option
view_option = 3;

%For 2D views, choose density or speed map underlay ('density', 'speed',
%'none')
map_underlay = 'none';

%Sun (the position of this one can be modified)
masses(1).name = 'Sun';
masses(1).mass = M_Sun;
masses(1).radii = 109.123/23455;
masses(1).x0 = 0;
masses(1).y0 = 0;
masses(1).z0 = 0;
masses(1).vx0 = 0;
masses(1).vy0 = 0;
masses(1).vz0 = 0;
masses(1).marker_RGB = [1,1,0];
masses(1).marker_size = 10;
masses(1).plot_trail = 0;
masses(1).init_orbit_M = NaN;
masses(1).init_orbit_P = 0;
masses(1).init_orbit_G = G;
masses(1).init_orbit_ecc = 0;
masses(1).init_orbit_s = 0;
masses(1).init_orbit_h0 = [0,0,0];
masses(1).init_orbit_hdot = [0,0,0];
masses(1).init_orbit_d = [1,0,0];
masses(1).init_orbit_alpha = 0;
masses(1).init_orbit_theta0 = 0;
masses(1).init_orbit_thetadot0 = 0;
masses(1).init_orbit_J = 0;
masses(1).init_orbit_j = 0;
masses(1).init_orbit_E = 0;

%Define names
names = {
    'Mercury',...
    'Venus',...
    'Earth',...
    'Mars',...
    'Jupiter',...
    'Saturn',...
    'Uranus',...
    'Neptune',...
    'Pluto'
    };

%Define masses of planets in Earth masses
M = [
    0.0553,...  %Mercury
    0.815,...   %Venus
    1.000,...   %Earth
    0.107,...   %Mars
    317.85,...  %Jupiter
    95.159,...  %Saturn
    14.500,...  %Uranus
    17.204,...  %Neptune
    0.003       %Pluto
    ];
%Convert to solar masses
M = M/332948;

%Semi-major axis of elliptical orbit in AU
a = [
    0.387,...   %Mercury
    0.723,...   %Venus
    1.000,...   %Earth
    1.523,...   %Mars
    5.202,...   %Jupiter
    9.575,...   %Saturn
    19.293,...  %Uranus
    30.246,...  %Neptune
    39.509      %Pluto
    ];

%Orbital eccentricity
ecc = [
    0.21,...   %Mercury
    0.01,...   %Venus
    0.02,...   %Earth
    0.09,...   %Mars
    0.05,...   %Jupiter
    0.06,...   %Saturn
    0.05,...   %Uranus
    0.01,...   %Neptune
    0.25       %Pluto
    ];

%Compute maximum separation from the sun
S = ( ecc + 1 ) .*a ;

%Orbital plane inclination (defined in degrees)
beta = [
    7.00,...   %Mercury
    3.39,...   %Venus
    0.00,...   %Earth
    1.85,...   %Mars
    1.31,...   %Jupiter
    2.49,...   %Saturn
    0.77,...   %Uranus
    1.77,...   %Neptune
    17.5       %Pluto
    ];
beta = beta*pi/180;

%Orbital ellipse semi-major axis direction
d = cell(1,length(M));
for n = 1:length(M)
    d{n} = [ cos(beta(n)), 0, sin(beta(n)) ];
end

%Rotation of ellipse about pointing direction
alpha = [
    0.0,...   %Mercury
    0.0,...   %Venus
    0.0,...   %Earth
    0.0,...   %Mars
    0.0,...   %Jupiter
    0.0,...   %Saturn
    0.0,...   %Uranus
    0.0,...   %Neptune
    0.0       %Pluto
    ];

%Object radii /Earth radii. Note 1 earth radii = 1/23,455 AU.
R = [
    0.383,...   %Mercury
    0.949,...   %Venus
    1.00,...    %Earth
    0.533,...   %Mars
    11.209,...  %Jupiter
    9.449,...   %Saturn
    4.007,...   %Uranus
    3.883,...   %Neptune
    0.187       %Pluto
    ];

%Plot trails
trails = [
    1,...   %Mercury
    1,...   %Venus
    1,...   %Earth
    1,...   %Mars
    1,...   %Jupiter
    1,...   %Saturn
    1,...   %Uranus
    1,...   %Neptune
    1       %Pluto
    ];

%RGB colour of orbits
C = { 
    (1/255)*[60,0,0],...         %Mercury
    (1/255)*[234,85,211],...     %Venus
    (1/255)*[106,175,255],...    %Earth
    (1/255)*[255,0,0],...        %Mars
    (1/255)*[255,128,100],...    %Jupiter
    (1/255)*[238,213,115],...    %Saturn
    (1/255)*[0,183,0],...        %Uranus
    (1/255)*[0,0,255],...        %Neptune
    (1/255)*[150,96,63]          %Pluto
    };

%Orbital initial phase angle /radians
theta0 = 2*pi*rand(1,length(M));

%Define a vector of planet masses
for n = 1:length(M)
    
    %Centre of mass initial position and velocity
    h0 = [0,0,0];
    hdot = [0,0,0];
    
    %Binary star circular orbit conditions
    [ux,vx1,vx2,...
        uy,vy1,vy2,...
        uz,vz1,vz2,...
        xc,x1,x2,...
        yc,y1,y2,...
        zc,z1,z2,...
        theta,theta_dot,t,j,J,E,P] =...
        two_body_elliptical_orbit(  G, M_Sun, M(n), S(n), ecc(n), theta0(n), h0, hdot, d{n}, alpha(n), 1 );
    
    %Sun (the position of this one can be modified)
    masses(n+1).name = names{n};
    masses(n+1).mass = M(n);
    masses(n+1).radii = R(n)/23455;
    masses(n+1).x0 = x2;
    masses(n+1).y0 = y2;
    masses(n+1).z0 = z2;
    masses(n+1).vx0 = vx2;
    masses(n+1).vy0 = vy2;
    masses(n+1).vz0 = vz2;
    masses(n+1).marker_RGB = C{n};
    masses(n+1).marker_size = 5;
    masses(n+1).plot_trail = trails(n);
    masses(n+1).init_orbit_M = M_Sun;
    masses(n+1).init_orbit_P = P;
    masses(n+1).init_orbit_G = G;
    masses(n+1).init_orbit_ecc = ecc(n);
    masses(n+1).init_orbit_s = S(n);
    masses(n+1).init_orbit_h0 = h0;
    masses(n+1).init_orbit_hdot = hdot;
    masses(n+1).init_orbit_d = d{n};
    masses(n+1).init_orbit_alpha = alpha(n);
    masses(n+1).init_orbit_theta0 = theta0(n);
    masses(n+1).init_orbit_thetadot0 = theta_dot;
    masses(n+1).init_orbit_J = J;
    masses(n+1).init_orbit_j = j;
    masses(n+1).init_orbit_E = E;
end

%Rings
rings = [];

%Clusters
clusters = [];

%End of code