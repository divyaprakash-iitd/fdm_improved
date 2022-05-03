function [solid,liquid,boundary,inner] = generate_flags(xp,yp,x,y)
    
    % Define flags
    nsolid = size(xp,1);
    solid = zeros(size(x,1),size(x,2),nsolid);
    for isolid = 1:nsolid
        [in, on] = inpolygon(x,y,xp(isolid,:),yp(isolid,:));
        solid(:,:,isolid) = logical(in + on);
    end

    liquid = ~sum(solid,3);

    boundary = sum(solid,3);

    % Rectangular domain
    boundary(1,:)   = 1;
    boundary(end,:) = 1;
    boundary(:,1)   = 1;
    boundary(:,end) = 1;

    boundary = logical(boundary);
    inner = ~boundary;
    solid = logical(solid);
end