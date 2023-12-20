function radius = findRadius(X,Y,Z)
    if length(X) == 0
        radius = 0;
    else
    arr = [X(:), Y(:), Z(:)];
    %D = pdist2(arr, arr);
    %radius = min(max(D, [], 2));
    meanX = mean(X(:));
    meanY = mean(Y(:));
    meanZ = mean(Z(:));
    centerPoint = [meanX, meanY, meanZ];
    distancesToCenter = sqrt(sum((arr - centerPoint).^2, 2));
    [sortedValues, indices] = sort(distancesToCenter, 'descend');
    %top100Values = sortedValues(1:100);
    radius = sortedValues(floor(length(X)*.01 + 1));
    %radius = max(distancesToCenter);
    end





end