function [solid,liquid,boundary,inner] = generate_flags(xp,yp,x,y)
    
    [in, on] = inpolygon(x,y,xp,yp);

    % Define flags
    solid = logical(in + on);
    liquid = ~solid;

    boundary = solid;
    boundary(1,:) = 1;
    boundary(end,:) = 1;
    boundary(:,1) = 1;
    boundary(:,end) = 1;
    boundary = logical(boundary);
    inner = ~boundary;
end