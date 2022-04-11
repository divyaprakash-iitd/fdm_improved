clear; clc; close all;

% calculate_coefficient()
% function A = calculate_coefficient()
    % Generate a rectangular grid

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
    t = linspace(0,2*pi,10*nx);
    cx = xmax/2; cy = ymax/2;
    radius = 0.2;
%     xp = radius*cos(t) + cx;
%     yp = radius*sin(t) + cy;
global a b
a = 0.2;
b = 0.3;
xp = a*cos(t) + cx;
yp = b*sin(t) + cy;
    
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

    % Intersection points
    xi = [];
    yi = [];
    % Loop over non-boundary nodes
    for i = 2:nx-1
        for j = 2:ny-1
            % Loop over only liquid nodes
            if liquid(i,j)
                % Check for all the four neighbouring nodes
                if solid(i+1,j)
                    % Calculate the intersection point
                    p = polyfit([x(i,j),x(i+1,j)],[y(i,j),y(i,j)],1);

                    t = atan2(y(i,j)-cy,x(i,j)-cx);
                    [xi, yi] = calculate_intersection_point(p,t,xi,yi,radius,cx,cy);
%                     plot(xi,yi,'mo')
                    x(i,j)
                   

                elseif solid(i-1,j)
                    p = polyfit([x(i,j),x(i-1,j)],[y(i,j),y(i,j)],1);
                    t = atan2(y(i,j)-cy,x(i,j)-cx);
                    [xi, yi] = calculate_intersection_point(p,t,xi,yi,radius,cx,cy);
%                     xqq = linspace(x(i,j),x(i-1,j),10);
%                     plot(xqq,polyval(p,xqq),'k-')
%                     plot(xi,yi,'mo')

                elseif solid(i,j+1)
                    p = polyfit([x(i,j),x(i,j)],[y(i,j),y(i,j+1)],1);
                    t = atan2(y(i,j)-cy,x(i,j)-cx);
                    [xi, yi] = calculate_intersection_point(p,t,xi,yi,radius,cx,cy);
%                     plot(xi,yi,'mo')

                elseif solid(i,j-1)
                    p = polyfit([x(i,j),x(i,j)],[y(i,j),y(i,j-1)],1);
                    t = atan2(y(i,j)-cy,x(i,j)-cx);
                    [xi, yi] = calculate_intersection_point(p,t,xi,yi,radius,cx,cy);
%                     plot(xi,yi,'mo')

                    
                end
            end
        end
    end
% end
plot(xi,yi,'mo')

function [xi, yi] = calculate_intersection_point(p,xo,xi,yi,r,cx,cy)
global a b
%     fun = @(x)  r*sin(x) + cy - p(1) * (r*cos(x) + cx) - p(2);
      fun = @(x)  b*sin(x) + cy - p(1) * (a*cos(x) + cx) - p(2);
    t = fsolve(fun,xo);
                  
%                     xi = [xi, r*cos(t)+cx];
%                     yi = [yi, r*sin(t)+cy];
                    xi = [xi, a*cos(t)+cx];
                    yi = [yi, b*sin(t)+cy];
end
