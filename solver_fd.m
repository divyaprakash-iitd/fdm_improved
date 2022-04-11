clear; clc; close all;

% Generate a rectangular grid

xmax = 2;
ymax = 1;
ar = xmax/ymax;

nx = 50;
ny = nx/ar;

% Create grid points
X = linspace(0,xmax,nx);
Y = linspace(0,ymax,ny);
[x,y] = ndgrid(X,Y);

% Generate a parametric closed curve
t = linspace(0,2*pi,nx);
cx = xmax/2; cy = ymax/2;
radius = 0.2;
xp = radius*cos(t) + cx;
yp = radius*sin(t) + cy;

[in, on] = inpolygon(x,y,xp,yp);

% Define flags
solid = logical(in + on);
liquid = ~solid;

hold on
plot(x(liquid),y(liquid),'rx')
plot(x(solid),y(solid),'bo')
plot(xp,yp,'-.')
axis off equal
