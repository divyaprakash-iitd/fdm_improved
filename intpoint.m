clear; clc; close all;

% Constants
LEFT        = 1;
RIGHT       = 2;
TOP         = 3;
BOTTOM      = 4;

xmax = 10;
ymax = 10;
ar = xmax/ymax;

nx = 200;
ny = nx/ar;

% Create grid points
X = linspace(0,xmax,nx);
Y = linspace(0,ymax,ny);
[x,y] = ndgrid(X,Y);

np = 60;
cx = xmax/2; cy = ymax/2;
[xp,yp] = closed_curve('cardoid',np,2,[cx,cy]);
% [xp,yp] = closed_curve('hypocycloid',np,3,[cx,cy]);
% [xp,yp] = closed_curve('tear',np,3,[cx,cy]);
% [xp,yp] = closed_curve('ellipse',np,3,[cx,cy]);
[xip,yip,xg,yg,beta] = calculate_coefficients(xp,yp,x,y);

[solid,liquid,boundary,inner] = generate_flags(xp,yp,x,y);  

T = zeros(size(x));

% Generate boundary conditions value

T(boundary) = 200;

Tb = abs(500*sin(x(solid)));
T(solid) = Tb;

for iter = 1:1000
    for i = 1:nx
        for j = 1:ny
            if inner(i,j)
                CC = 2*(1/beta(i,j,LEFT)/beta(i,j,RIGHT) + 1/beta(i,j,TOP)/beta(i,j,BOTTOM));
                CL = 2/beta(i,j,LEFT)*(1/(beta(i,j,LEFT) + beta(i,j,RIGHT)));
                CR = 2/beta(i,j,RIGHT)*(1/(beta(i,j,LEFT) + beta(i,j,RIGHT)));
                CT = 2/beta(i,j,TOP)*(1/(beta(i,j,TOP) + beta(i,j,BOTTOM)));
                CB = 2/beta(i,j,BOTTOM)*(1/(beta(i,j,TOP) + beta(i,j,BOTTOM)));
    
                T(i,j) = 1/CC * (T(i+1,j)*CR + T(i-1,j)*CL + T(i,j-1)*CB + T(i,j+1)*CT);
    
            end
        end
    end
end

contourf(xg,yg,T,50,'edgecolor','none')
colormap("jet")
axis equal
colorbar
hold on
plot(xp,yp,'LineWidth',2)
plot(xip,yip,'ko','MarkerSize',6)
mesh(x,y,ones(size(x)),'facealpha',0,'edgealpha',0.2)

% figure()
% scatter(xg(:),yg(:),T(:),T(:),'.','LineWidth',4)

% plot(x    g,yg,'w+')
% hold on
% plot(x(liquid),y(liquid),'rx')
% plot(x(solid),y(solid),'bx')
% plot(xp,yp,'k-.')
% axis equal
% grid on
% plot(xip,yip,'r*')