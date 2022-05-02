function [xp,yp] = closed_curve(name,np,l,center)
    % Generate a parametric closed curve
    t = linspace(0,2*pi,np);
    cx = center(1); cy = center(2);
    if strcmp(name,'ellipse')
        a = l;
        b = a/ar;
        xp = a* cos(t);
        yp = b*sin(t);
    end

    if strcmp(name,'hypocycloid') 
        xp = l*cos(t)+cos(l*t);
        yp = l*sin(t)-sin(l*t);
    end

    xp = xp + cx;
    yp = yp + cy;
end