%vrot
% Three dimensional rotation function, using a generalized version of
% Rodrigues' formula.
%
% LAST UPDATED by Andy French Jan 2013
%
% Syntax: [x,y,z] = vrot(x,y,z,wx,wy,wz,ax,ay,az,theta)
%
% x,y,z are NxM matrices of Cartesian (x,y,z) coordinates to be rotated
% [wx,wy,wz] is the Cartesian vector describing the axis of rotation
% [ax,ay,az] is the Cartesian vector describing the centre of the rotation
% theta is the anti-clockwise angle of rotation in radians about the rotation axis

function [x,y,z] = vrot(x,y,z,wx,wy,wz,ax,ay,az,theta)

%Determine length of rotation axis vector
W = sqrt( wx^2 + wy^2 + wz^2 );

%Get dimensions of input matrices
[N,M] = size(x);

%Construct 3D matrices of coordinates. The third (pages) dimension
%corresponds to the Cartesian x,y,z spatial dimensions
r = zeros(N,M,3);
a = ones(N,M,3);
w = ones(N,M,3);
r(:,:,1) = x;
r(:,:,2) = y;
r(:,:,3) = z;
a(:,:,1) = ax*a(:,:,1);
a(:,:,2) = ay*a(:,:,1);
a(:,:,3) = az*a(:,:,1);
w(:,:,1) = wx*w(:,:,1);
w(:,:,2) = wy*w(:,:,2);
w(:,:,3) = wz*w(:,:,3);

%Compute new coordinates via the generalized Rodrigues formula
new_r = a*(1-cos(theta)) +...
    r*cos(theta) +...
    cross( w, r-a, 3 ) * sin(theta) / W +...
    ( 1-cos(theta) )*repmat( ( dot(r,w,3) - dot(a,w,3) ),[1,1,3] ).*w / (W^2) ;
    
%Extract outputs from new_r
x(:,:) = new_r(:,:,1);
y(:,:) = new_r(:,:,2);
z(:,:) = new_r(:,:,3);

%End of code