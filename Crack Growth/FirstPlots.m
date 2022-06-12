% Plot symbols at all the input points
for i = 1:Timesteps
    Pointnumber = PointsSize(i);
    k = mod(i-1,6)+1;
    for j = 1:Pointnumber
        testpoint = Points{i,j};
        plot3(testpoint(1), testpoint(2), testpoint(3),symbolchange{k},'MarkerEdgeColor','k',...
            'MarkerFaceColor',colourchange{k},'MarkerSize',12.5);
    end
end

%Plot lines connecting input points, and calculate average point spacing
distance = 0;
for i = 1:Timesteps
    Pointnumber = PointsSize(i);
    k = mod(i-1,6)+1;
    linedistance = 0;
    for j = 1:Pointnumber-1
        testpoint = Points{i,j};
        testpoint2 = Points{i,j+1};
        linedistance = linedistance + norm(testpoint2-testpoint);
        plot3([testpoint(1),testpoint2(1)],[testpoint(2),testpoint2(2)],[testpoint(3),...
            testpoint2(3)],'color',colourchange{k}.*0.7);
    end
    distance = distance + linedistance/(Pointnumber-1);
end
distance = distance/Timesteps;
disp(sprintf('Average point spacing %8.5E', distance ));