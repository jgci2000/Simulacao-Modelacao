%Gravity sim default inputs - with cluster!

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

%Set dimensions of space
Lmax = 8;

%Timestep
dt = 0.1;

%Gravitation exponent
Ng = 2;

%Repulsion exponent
Nr = 5;

%Repulsion magnitude
RM = 0.01;

%Star masses in Solar masses
M1 = 1;
M2 = 1;

%Initial star separation /AU
s = 6;

%Initial phase (polar angle) of two-mass elliptical orbit
theta0 = 0;

%Elliptical orbit eccentricity
ecc = 0;

%Orientation of semi-major axis of ellipse
d = [1,0,0];

%Rotation of ellipse clockwise about d /radians
alpha = 0;

%Centre of mass initial coordinates
h0 = [0,0,0];

%Cetre of mass velocity
hdot = [0,0,0];

%Boundary coefficient of restitution (if [] then there is
%no boundary)
k = 0.5;

%Initial view option
view_option = 2;

%For 2D views, choose density or speed map underlay ('density', 'speed',
%'none')
map_underlay = 'none';

%Binary star circular orbit conditions
[ux,vx1,vx2,...
    uy,vy1,vy2,...
    uz,vz1,vz2,...
    xc,x1,x2,...
    yc,y1,y2,...
    zc,z1,z2,...
    theta,theta_dot,t,j,J,E,P] =...
    two_body_elliptical_orbit(  G, M1, M2, s, ecc, theta0, h0, hdot, d, alpha, 1 );

%Star #1 (This can be modified dynamically via the user)
masses(1).name = 'Sun #1';
masses(1).mass = M1;
masses(1).radii = 0.1;
masses(1).x0 = x1;
masses(1).y0 = y1;
masses(1).z0 = z1;
masses(1).vx0 = vx1;
masses(1).vy0 = vy1;
masses(1).vz0 = vz1;
masses(1).marker_RGB = [];
masses(1).marker_size = 5;
masses(1).plot_trail = 1;
masses(1).init_orbit_M = M2;
masses(1).init_orbit_P = P;
masses(1).init_orbit_G = G;
masses(1).init_orbit_ecc = ecc;
masses(1).init_orbit_s = s;
masses(1).init_orbit_h0 = h0;
masses(1).init_orbit_hdot = hdot;
masses(1).init_orbit_d = d;
masses(1).init_orbit_alpha = alpha;
masses(1).init_orbit_theta0 = theta0;
masses(1).init_orbit_thetadot0 = theta_dot;
masses(1).init_orbit_J = J;
masses(1).init_orbit_j = j;
masses(1).init_orbit_E = E;

%Star #2
masses(2).name = 'Sun #2';
masses(2).mass = M2;
masses(2).radii = 0.1;
masses(2).x0 = x2;
masses(2).y0 = y2;
masses(2).z0 = z2;
masses(2).vx0 = vx2;
masses(2).vy0 = vy2;
masses(2).vz0 = vz2;
masses(2).marker_RGB = [];
masses(2).marker_size = 5;
masses(2).plot_trail = 1;
masses(2).init_orbit_M = M1;
masses(2).init_orbit_P = P;
masses(2).init_orbit_G = G;
masses(2).init_orbit_ecc = ecc;
masses(2).init_orbit_s = s;
masses(2).init_orbit_h0 = h0;
masses(2).init_orbit_hdot = hdot;
masses(2).init_orbit_d = d;
masses(2).init_orbit_alpha = alpha;
masses(2).init_orbit_theta0 = theta0+pi;
masses(2).init_orbit_thetadot0 = theta_dot;
masses(2).init_orbit_J = J;
masses(2).init_orbit_j = j;
masses(2).init_orbit_E = E;

%Rings (1)
rings(1).xc = x1;
rings(1).yc = y1;
rings(1).zc = z1;
rings(1).vxc = vx1;
rings(1).vyc = vy1;
rings(1).vzc = vz1;
rings(1).num_rings = 30;
rings(1).arc_separation_AU = 1*pi/30;
rings(1).first_ring_radius_AU = 2;
rings(1).ring_radius_diff_AU = 0.1;
rings(1).d = [1,0,0];
rings(1).alpha = 0;
rings(1).mass_at_centre = M1;
rings(1).marker_RGB = [0,0,0];

%Rings (2)
rings(2).xc = x2;
rings(2).yc = y2;
rings(2).zc = z2;
rings(2).vxc = vx2;
rings(2).vyc = vy2;
rings(2).vzc = vz2;
rings(2).num_rings = 30;
rings(2).arc_separation_AU = 1*pi/30;
rings(2).first_ring_radius_AU = 2;
rings(2).ring_radius_diff_AU = 0.1;
rings(2).d = [1,0,0];
rings(2).alpha = 0;
rings(2).mass_at_centre = M2;
rings(2).marker_RGB = [1,0,1];

%Clusters
clusters = [];

%End of code