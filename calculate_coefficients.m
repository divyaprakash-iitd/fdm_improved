clear; clc; close all;

xmax = 1;
ymax = 1;
ar = xmax/ymax;

nx = 20;
ny = nx/ar;

% Create grid points
X = linspace(0,xmax,nx);
Y = linspace(0,ymax,ny);
[x,y] = ndgrid(X,Y);

% Generate a parametric closed curve
np = 60;
t = linspace(0,2*pi,np);
cx = xmax/2; cy = ymax/2;
radius = 0.2;

a = 0.25;
b = 0.25;
xp = a* cos(t) + cx;
yp = b*sin(t) + cy;

% Tells us the index of x and y which are inside or on the closed curve.so                          
[in, on] = inpolygon(x,y,xp,yp);

% Define flags
solid = logical(in + on);
liquid = ~solid;

hold on
plot(x(liquid),y(liquid),'rx')
plot(x(solid),y(solid),'bx')
plot(xp,yp,'k-.')
axis equal
grid on

% Loop over non-boundary nodes
for i = 2:nx-1
    for j = 2:ny-1
        % Loop over only liquid nodes
        if liquid(i,j)
%             plot(x(i,j),y(i,j),'ro')
            % Check for all the four neighbouring nodes
            if solid(i+1,j)
%                 plot(x(i+1,j),y(i+1,j),'ro')
                % Calculate the intersection point
                px = [x(i,j),x(i+1,j)];
                py = [y(i,j),y(i+1,j)];
                [xi, yi] = polyxpoly(xp,yp,px,py);     
                plot(xi,yi,'mo')
            end
                
            if solid(i-1,j)
%                 plot(x(i-1,j),y(i-1,j),'ro')
                px = [x(i-1,j),x(i,j)];
                py = [y(i-1,j),y(i,j)];
                [xi, yi] = polyxpoly(xp,yp,px,py);
                plot(xi,yi,'mo')
            end
                
            if solid(i,j+1)
%                 plot(x(i,j+1),y(i,j+1),'ro')
                px = [x(i,j),x(i,j+1)];
                py = [y(i,j),y(i,j+1)];
                [xi, yi] = polyxpoly(xp,yp,px,py);
                plot(xi,yi,'mo')
            end
                           
            if solid(i,j-1)
%                 plot(x(i,j-1),y(i,j-1),'ro')
                px = [x(i,j-1),x(i,j)];
                py = [y(i,j-1),y(i,j)];
                [xi, yi] = polyxpoly(xp,yp,px,py);
                plot(xi,yi,'mo')             
            end
        end
    end
end
