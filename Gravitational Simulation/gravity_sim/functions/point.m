%point
% Rotate a matrix of position vectors such that a particular position
% vector rotates to a specified direction. For example, the coordinates of
% the vertices of an aircraft model which want to be rotated such that the
% aircraft faces a given direction without unintended roll about that axis.
%
% LAST UPDATED by Andy French Jan 2013
%
% Syntax: [x,y,z] = point(x,y,z,bx,by,bz,wx,wy,wz)
%
% x,y,z are NxM matrices of Cartesian (x,y,z) coordinates to be rotated.
% [bx,by,bz] is the Cartesian vector describing vector in the scene which
% needs to be rotated to [wx,wy,wz].
%
% Note all rotations are about the origin, e.g. the centre of mass of the
% object being rotated. If the origin is currently at an undesirable place,
% shift it before using this function.

function [x,y,z] = point(x,y,z,bx,by,bz,wx,wy,wz)

%Determine length of input vectors
W = sqrt( wx^2 + wy^2 + wz^2 );
B = sqrt( bx^2 + by^2 + bz^2 );

%Form 'average' rotation axis
ww = [wx,wy,wz]/W + [bx,by,bz]/B;
WW = sqrt( ww(1)^2 + ww(2)^2 + ww(3)^2 );
wwx = ww(1);
wwy = ww(2);
wwz = ww(3);

%Get dimensions of input matrices
[N,M] = size(x);

%Construct 3D matrices of coordinates. The third (pages) dimension
%corresponds to the Cartesian x,y,z spatial dimensions
r = zeros(N,M,3);
a = zeros(N,M,3);  %Note all rotations are about the origin. 
w = ones(N,M,3);
r(:,:,1) = x;
r(:,:,2) = y;
r(:,:,3) = z;
w(:,:,1) = wwx*w(:,:,1);
w(:,:,2) = wwy*w(:,:,2);
w(:,:,3) = wwz*w(:,:,3);

%Compute new coordinates via the generalized Rodrigues formula
theta = pi;
r = a*(1-cos(theta)) +...
    r*cos(theta) +...
    cross( w, r-a, 3 ) * sin(theta) / WW +...
    ( 1-cos(theta) )*repmat( ( dot(r,w,3) - dot(a,w,3) ),[1,1,3] ).*w / (WW^2) ;

%Now rotate by pi about [wx,wy,wz] to correct for the problem that the new
%scene is upside down
w = ones(N,M,3);
w(:,:,1) = wx*w(:,:,1);
w(:,:,2) = wy*w(:,:,2);
w(:,:,3) = wz*w(:,:,3);
theta = pi;
r = a*(1-cos(theta)) +...
    r*cos(theta) +...
    cross( w, r-a, 3 ) * sin(theta) / W +...
    ( 1-cos(theta) )*repmat( ( dot(r,w,3) - dot(a,w,3) ),[1,1,3] ).*w / (W^2) ;
    
%Extract outputs from r
x(:,:) = r(:,:,1);
y(:,:) = r(:,:,2);
z(:,:) = r(:,:,3);

%End of code