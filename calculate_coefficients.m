function [xip,yip,xg,yg,beta] = calculate_coefficients(xp,yp,x,y)   

                 
[in, on] = inpolygon(x,y,xp,yp);

% Define flags
solid = logical(in + on);
liquid = ~solid;

% hold on
% plot(x(liquid),y(liquid),'rx')
% plot(x(solid),y(solid),'bx')
% plot(xp,yp,'k-.')
% axis equal
% grid on

xip = [];
yip = [];

nx = size(x,1);
ny = size(x,2);

xg = x; yg = y;
xing(solid) = nan;
yg(solid) = nan;

% Coefficient matrix
% Constants
LEFT        = 1;
RIGHT       = 2;
TOP         = 3;
BOTTOM      = 4;

beta = ones(nx,ny,4);

% Loop over non-boundary nodes
for i = 2:nx-1
    for j = 2:ny-1
        % Loop over only liquid nodes
        if liquid(i,j)
            if solid(i+1,j)
                px = [x(i,j),x(i+1,j)];
                py = [y(i,j),y(i+1,j)];
                [xi, yi] = polyxpoly(xp,yp,px,py);     
%                 plot(xi,yi,'mo')
                xip = [xip, xi];
                yip = [yip, yi];

                xg(i+1,j) = xi;
                yg(i+1,j) = yi;

                % Right
                % Distance between intersection point and the grid point
                beta(i,j,RIGHT) = norm([xi-x(i+1,j),yi-y(i+1,j)]);

            end
                
            if solid(i-1,j)
                px = [x(i-1,j),x(i,j)];
                py = [y(i-1,j),y(i,j)];
                [xi, yi] = polyxpoly(xp,yp,px,py);
%                 plot(xi,yi,'mo')
                xip = [xip, xi];
                yip = [yip, yi];

                xg(i-1,j) = xi;
                yg(i-1,j) = yi;

                beta(i,j,LEFT) = norm([xi-x(i-1,j),yi-y(i-1,j)]);
            end
                
            if solid(i,j+1)
                px = [x(i,j),x(i,j+1)];
                py = [y(i,j),y(i,j+1)];
                [xi, yi] = polyxpoly(xp,yp,px,py);
%                 plot(xi,yi,'mo')
                xip = [xip, xi];
                yip = [yip, yi];

                xg(i,j+1) = xi;
                yg(i,j+1) = yi;

                beta(i,j,TOP) = norm([xi-x(i,j+1),yi-y(i,j+1)]);
            end
                           
            if solid(i,j-1)
                px = [x(i,j-1),x(i,j)];
                py = [y(i,j-1),y(i,j)];
                [xi, yi] = polyxpoly(xp,yp,px,py);
                xip = [xip, xi];
                yip = [yip, yi];
%                 plot(xi,yi,'mo') 
                
                xg(i,j-1) = xi;
                yg(i,j-1) = yi;

                beta(i,j,BOTTOM) = norm([xi-x(i,j-1),yi-y(i,j-1)]);
            end
        end
    end
end
end
