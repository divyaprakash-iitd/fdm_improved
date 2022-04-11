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
[x,y] = meshgrid(X,Y);

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
plot(xp,yp,'k-.')
axis off equal

% Solve a dirichlet heat conduction problem
% Boundary condition
Tw = 10; % West
Te = 20; % East
Tn = 30;  % North
Ts = 40; % South


% Description: Solves a 2D steady state heat transfer problem with no heat
% generation, with Dirichlet boundary conditions at the boundaries.

% Constants
N           = 5;
LEFT        = 1;
RIGHT       = 2;
TOP         = 3;
BOTTOM      = 4;
CENTER      = 5;

% Grid Size
dx  = xmax/(nx-1);
dy  = ymax/(ny-1);

Nx  = nx;
Ny  = ny;

% Coefficient Matrix
A = ones(Nx,Ny,N);
A(:,:,CENTER) = -4;

% RHS Matrix
B = zeros(Nx,Ny);

% Solution Matrix
Ti = 300; % Initial Guess
T = Ti*ones(Nx,Ny);


% Dirichlet BC (Left and Right)
T(:,1)      = Te;
T(:,end)    = Tw;
% Dirichlet BC (Top and Bottom)
T(1,:)      = Tn;
T(end,:)    = Ts;


T = solve_poisson(A,B,T);

% Implement a function which calculates the coefficients
A = calculate_coefficients(x,y,xp,yp);