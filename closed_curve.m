function [xp,yp] = closed_curve(name,np,l,center)
    % Generate a parametric closed curve
    t = linspace(0,2*pi,np);
    cx = center(1); cy = center(2);
    if strcmp(name,'ellipse')
        a = l;
        ar = 2;
        b = a/ar;
        xp = a* cos(t);
        yp = b*sin(t);
    end

    if strcmp(name,'hypocycloid') 
        xp = l*cos(t)+cos(l*t);
        yp = l*sin(t)-sin(l*t);
    end

    if strcmp(name,'tear') 
        m = 1;
        xp = cos(t);
        yp = sin(t).*(sin(1/2*t)).^m;
    end

    if strcmp(name,'cardoid') 
        m = 1;
        xp = l*cos(t).*(1-cos(t));
        yp = l*sin(t).*(1-cos(t));
    end

    xp = xp + cx;
    yp = yp + cy;
end