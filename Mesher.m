% Body Parameters
dSample = 0.5;
lSample = 5.3;

% Tip Parameters
tipAngle = deg2rad(50.0);

% Shank Parameters
dShank = 1.25;
lShank = 0.4;

% Thread Parameters
threadCount = 4.0;
threadDepth = 0.125;
approachAngle = deg2rad(22.5);
recedingAngle = deg2rad(22.5);
threadCrest = 0.03;

% Mesh Parameters
nFace = 3;
nSlice = round((nFace*lSample)/(pi*dSample));

[verts1, verts2, verts3, centroids, threadHeightMap, tipLength] = createCylinderArray(nSlice, nFace, dSample, tipAngle,...
    lSample, dShank, lShank, threadCount, threadDepth, approachAngle, recedingAngle, threadCrest);
[nPolygons, polygonWidth, polygonArea] = getMeshParameters(dSample, nFace, nSlice, verts1, verts2, verts3);

movieFig = figure('Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
movieAxes = axes(movieFig);
movieAxes.XLim = [-lSample/1.75,lSample/1.75];
movieAxes.YLim = [-lSample/1.75,lSample/1.75];
movieAxes.ZLim = [-lSample/1.75,lSample/1.75];
movieAxes.DataAspectRatio = [1 1 1];
movieAxes.XLabel.String = 'X-Axis [cm]';
movieAxes.YLabel.String = 'Y-Axis [cm]';
movieAxes.ZLabel.String = 'Z-Axis [cm]';
movieAxes.FontSize = 16;
movieAxes.XGrid = 'on';
movieAxes.YGrid = 'on';
movieAxes.ZGrid = 'on';
movieAxes.View = [-45 20];
patch(movieAxes, [verts1(1,:); verts2(1,:); verts3(1,:)],...
    [verts1(2,:); verts2(2,:); verts3(2,:)],...
    [verts1(3,:); verts2(3,:); verts3(3,:)],...
    [0.70, 0.70, 0.70], 'EdgeAlpha', '0.1');

function [verts1, verts2, verts3, centroids, threadHeightMap, tipLength] = createCylinderArray(nSlice, nFace, dSample, tipAngle,...
    lSample, dShank, lShank, threadCount, threadDepth, approachAngle, recedingAngle, threadCrest)
   
    polygonWidth = (nFace * dSample * sin(pi / nFace)) / (nFace+1);
    [flatPoints, pointSampleIndex] = getFlatPoints(nSlice, nFace, polygonWidth, lSample);
    threadHeightMap = getThreadHeightMap(threadCount, threadDepth, approachAngle, recedingAngle, threadCrest, dSample, lSample, pointSampleIndex, flatPoints);

    %Initialize mesh variables
    X = (zeros(3,2*(nSlice-1)*nFace));
    Y = (zeros(3,2*(nSlice-1)*nFace));
    Z = (zeros(3,2*(nSlice-1)*nFace));
    points = getPoints(nSlice, nFace, dSample, tipAngle, lSample, dShank, lShank, threadHeightMap);

    currPoint = 1;
    for j = 1:(nSlice-1)
        for k = 1:2*nFace
            spacer = (j-1)*(nFace+1);
            if k <= nFace
                p1 = points(:, k + spacer);
                p2 = points(:, k + 1 + spacer);
                p3 = points(:, k + nFace + 2 + spacer);
            else
                p1 = points(:, (k - nFace) + spacer);
                p2 = points(:, k + 1 + spacer);
                p3 = points(:, k + 2 + spacer);
            end
            X(:,currPoint) = [p1(1,1);p2(1,1);p3(1,1)];
            Y(:,currPoint) = [p1(2,1);p2(2,1);p3(2,1)];
            Z(:,currPoint) = -[p1(3,1);p2(3,1);p3(3,1)];
            currPoint = currPoint + 1;
        end
    end

    if tipAngle == pi
        tipLength = 0;
    else
        tipLength = dSample / (2*tan(tipAngle/2));
    end
    verts1 = [X(1,:); Y(1,:); Z(1,:)];
    verts2 = [X(2,:); Y(2,:); Z(2,:)];
    verts3 = [X(3,:); Y(3,:); Z(3,:)];
    centroids = [mean(X); mean(Y); mean(Z)];
    
end

function [flatPoints, pointSampleIndex] = getFlatPoints(nSlice, nFace, polygonWidth, lSample)

    flatPoints = zeros(2,nSlice*(nFace+1));
    pointSampleIndex = zeros(1,nSlice*(nFace+1));

    totalWidth = polygonWidth*((nFace+1)-1);
    currPoint = 1;

    startX = -totalWidth/2;
    startY = -lSample/2;

    for currSlice = 1:nSlice
        currY = startY + (currSlice-1)*lSample/(nSlice-1);

        for currFace = 1:(nFace+1)
            currX = startX + (currFace-1)*polygonWidth;

            flatPoints(:,currPoint) = [currX; currY];
            pointSampleIndex(1,currPoint) = 1;
            currPoint = currPoint + 1;
        end
    end

end

function threadHeightMap = getThreadHeightMap(threadCount, threadDepth, approachAngle, recedingAngle, threadCrest, dSample, lSample, pointSampleIndex, flatPoints)

    %Initialize the parameters used to define the geometric positions of
    %the threads
    nThreads = floor(lSample * threadCount);
    threadSlope = 1 / (pi * dSample * threadCount);
    approachDiameter = threadDepth*tan(approachAngle);
    recedingDiameter = threadDepth*tan(recedingAngle);
    inclinedApproachDiameter = approachDiameter/cos(threadSlope);
    inclinedCrestDiameter = threadCrest/cos(threadSlope);
    inclinedRecedingDiameter = recedingDiameter/cos(threadSlope);
    totalInclinedDiameter = inclinedApproachDiameter + inclinedCrestDiameter + inclinedRecedingDiameter;

    %Solve the geometric thread positions
    threadIntercepts = zeros(nThreads, 4);
    for i = 1:nThreads
        midLine = (i - 1) * 1/threadCount - lSample/2;
        approachLine = midLine + totalInclinedDiameter/2;
        approachCrestLine = midLine + totalInclinedDiameter/2 - inclinedApproachDiameter;
        recedingCrestLine = midLine - totalInclinedDiameter/2 + inclinedRecedingDiameter;
        recedingLine = midLine - totalInclinedDiameter/2;
        threadIntercepts(i,:) = [approachLine, approachCrestLine, recedingCrestLine, recedingLine];
    end

    %Solve for the left edges (min X value) of the samples
    sampleFirstPoint = zeros(1,length(pointSampleIndex(1,:)));
    for i = 1:length(pointSampleIndex)-1
        if pointSampleIndex(i+1) ~= pointSampleIndex(i)
            sampleFirstPoint(i+1) = 1;
        end
    end
    sampleFirstPoint(1,1) = 1;
    xCoords = flatPoints(1,:);
    yCoords = flatPoints(2,:);
    sampleStart = xCoords(sampleFirstPoint(1,:)==1);

    %Determine whether each polygon falls within the thread. If they do,
    %apply the proper thread angle fluence correction factor
    threadHeightMap = (zeros(1,length(xCoords(1,:))));

    for currPoint = 1 : length(xCoords(1,:))
        yConst = yCoords(currPoint);
        xConst = threadSlope * (xCoords(currPoint) - sampleStart(pointSampleIndex(currPoint)));
        for currThread = 1 : nThreads

            eqn1 = (yConst) <=  xConst + threadIntercepts(currThread, 1);
            eqn2 = (yConst) >=  xConst + threadIntercepts(currThread, 2);
            eqn3 = (yConst) <=  xConst + threadIntercepts(currThread, 3);
            eqn4 = (yConst) >=  xConst + threadIntercepts(currThread, 4);
            eqn5 = (yConst) <=  xConst + threadIntercepts(currThread, 2);
            eqn6 = (yConst) >=  xConst + threadIntercepts(currThread, 3);

            if eqn1 && eqn2
                A = -threadSlope;
                B = 1;
                C = (threadSlope*sampleStart(pointSampleIndex(currPoint))-threadIntercepts(currThread, 1));
                x1 = xCoords(currPoint);
                y1 = yCoords(currPoint);
                distance = abs((A*x1+B*y1+C)/sqrt(A^2+B^2));
                threadHeightMap(1,currPoint) = distance / tan(approachAngle);
            elseif eqn3 && eqn4
                A = -threadSlope;
                B = 1;
                C = (threadSlope*sampleStart(pointSampleIndex(currPoint))-threadIntercepts(currThread, 4));
                x1 = xCoords(currPoint);
                y1 = yCoords(currPoint);
                distance = abs((A*x1+B*y1+C)/sqrt(A^2+B^2));
                threadHeightMap(1,currPoint) =  distance / tan(recedingAngle);
            elseif eqn5 && eqn6
                threadHeightMap(1,currPoint) =  threadDepth;
            end

        end
    end

end

function points = getPoints(nSlice, nFace, dSample, tipAngle, lSample, dShank, lShank, threadHeightMap)
    
    currPoint = 1;
    points = zeros(3,nSlice * (nFace+1));
    taperTerminationZ = lSample/2 - dSample/(2*tan(tipAngle/2));

    for currRow = 1:nSlice
        for currColumn =  1:nFace+1

            currZ = -1 * ((lSample/2) - (currRow-1) * (lSample/(nSlice-1)));

            %Taper
            if currZ >= taperTerminationZ && tipAngle ~= pi && currZ ~= lSample/2
                currX = ((lSample/2-currZ)/(lSample/2-taperTerminationZ))*(dSample/2 + threadHeightMap(currPoint))*cos(2*(currColumn-1)*pi/nFace);
                currY = ((lSample/2-currZ)/(lSample/2-taperTerminationZ))*(dSample/2 + threadHeightMap(currPoint))*sin(2*(currColumn-1)*pi/nFace);
            
            %Tip
            elseif currZ >= taperTerminationZ && tipAngle ~= pi && currZ == lSample/2
                tipZ = -1 * ((lSample/2) - (currRow-1.01) * (lSample/(nSlice-1)));
                currX = ((lSample/2-tipZ)/(lSample/2-taperTerminationZ))*(dSample/2 + threadHeightMap(currPoint))*cos(2*(currColumn-1)*pi/nFace);
                currY = ((lSample/2-tipZ)/(lSample/2-taperTerminationZ))*(dSample/2 + threadHeightMap(currPoint))*sin(2*(currColumn-1)*pi/nFace);
            
            %Shank
            elseif currZ < -lSample/2 + lShank
                currX = (dShank/2)*cos(2*(currColumn-1)*pi/nFace);
                currY = (dShank/2)*sin(2*(currColumn-1)*pi/nFace);
            
            %Thread body
            else
                currX = (dSample/2 + threadHeightMap(currPoint))*cos(2*(currColumn-1)*pi/nFace);
                currY = (dSample/2 + threadHeightMap(currPoint))*sin(2*(currColumn-1)*pi/nFace);
            end

            if abs(currX) <= 1.0e-6
                currX = 0.0;
            end
            if abs(currY) <= 1.0e-6
                currY = 0.0;
            end
            if abs(currZ) <= 1.0e-6
                currZ = 0.0;
            end
            points(:, currPoint) = [currX; currY; currZ];

            currPoint = currPoint + 1;
        end
    end
    
end

function [nPolygons, polygonWidth, polygonArea] = getMeshParameters(dSample, nFace, nSlice, verts1, verts2, verts3)
    
    nPolygons = 2 * nFace * (nSlice - 1);

    polygonArea = zeros(1,nPolygons);
    for currPolygon = 1:nPolygons
        polygonArea(1,currPolygon) = 0.5*norm(cross((verts2(:,currPolygon)-verts1(:,currPolygon)),(verts3(:,currPolygon)-verts1(:,currPolygon))));
    end

    polygonWidth = (nFace * dSample * sin(pi / nFace)) / (nFace+1);
    
end