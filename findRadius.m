function radius = findRadius(X,Y,Z)

    arr = [X(:), Y(:), Z(:)];
    %D = pdist2(arr, arr);
    %radius = min(max(D, [], 2));
    meanX = mean(X(:));
    meanY = mean(Y(:));
    meanZ = mean(Z(:));
    centerPoint = [meanX, meanY, meanZ];
    distancesToCenter = sqrt(sum((arr - centerPoint).^2, 2));
    radius = max(distancesToCenter);






end