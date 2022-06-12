%Resample long lines
for i = 1:Timesteps
    Pointnumber = PointsSize(i);
    newpoint = 0;
    n = mod(i-1,6)+1;
    for j = 1:Pointnumber-1
        testpoint = Points{i,j};
        testpoint2 = Points{i,j+1};
        localdistance = norm(testpoint2-testpoint);
        
        newpoint = newpoint + 1;
        InterPoints{i,newpoint} = Points{i,j};
        if localdistance > (distance/SamplingSteps)
            m = ceil(SamplingSteps*localdistance/distance);
            %Add more new points
            for k = [1:m-1]
                newpoint = newpoint + 1;
                InterPoints{i,newpoint} = testpoint*(1-k/m) + testpoint2*k/m;
                plotpoint =  InterPoints{i,newpoint};
                plot3(plotpoint(1), plotpoint(2), plotpoint(3),symbolchange{n},'MarkerEdgeColor','k',...
                'MarkerFaceColor',colourchange{n}/2,'MarkerSize',5);

            end
        end
    end
    newpoint = newpoint + 1;
    InterPoints{i,newpoint} = Points{i,Pointnumber};
    InterPointsSize(i) = newpoint;
end