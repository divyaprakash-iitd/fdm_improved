clear; clc; close all;

% calculate_coefficient()
% function A = calculate_coefficient()
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
    plot(x(liquid),y(liquid),'r.')
    plot(x(solid),y(solid),'b.')
    plot(xp,yp,'k-.')
    axis equal
    grid on

    opts = optimoptions(@fsolve,'Algorithm', 'levenberg-marquardt');
    xx = oncurve(xp,yp,xp(1),yp(1));
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
                    [xi, yi] = calculate_intersection_point(p,10,xi,yi,radius,cx,cy);
                    xqq = linspace(x(i,j),x(i+1,j),10);
                    plot(xqq,polyval(p,xqq),'-')
                    plot(xi,yi,'kx')
                   
%                     return

                elseif solid(i-1,j)
                    p = polyfit([x(i,j),x(i-1,j)],[y(i,j),y(i,j)],1);
                    [xi, yi] = calculate_intersection_point(p,x(i,j),xi,yi);
                    [xi, yi] = calculate_intersection_point(p,10,xi,yi,radius,cx,cy);
                    xqq = linspace(x(i,j),x(i+1,j),10);
                    plot(xqq,polyval(p,xqq),'-')
                    plot(xi,yi,'kx')
% 
%                 elseif solid(i,j+1)
%                      p = polyfit([x(i,j),x(i,j)],[y(i,j),y(i,j+1)],1);
%                      [xi, yi] = calculate_intersection_point(p,y(i,j),xi,yi);
% 
%                 elseif solid(i,j-1)
%                     p = polyfit([x(i,j),x(i,j)],[y(i,j),y(i,j-1)],1);
%                     [xi, yi] = calculate_intersection_point(p,y(i,j),xi,yi);

                end
            end
        end
    end
% end
function [xi, yi] = calculate_intersection_point(p,xo,xi,yi,r,cx,cy)
    t = fzero(@(x) r*sin(x) + cy - p(1) * (r*cos(x) + cx) - p(2),xo);
                    xi = [xi, r*cos(t)+cx];
                    yi = [yi, r*sin(t)+cy];
end

function xx = oncurve(xp,yp,x,y)
    [~, on] = inpolygon(xp,yp,x,y);
    if any(on)
        xx = any(on) - 1;
    else
        xx = ~any(on);
    end
    xx = double(xx);
end