handlevelocities = figure;
hold on
daspect([1 1 1])

for i = 1:Timesteps-1
    Pointnumber = PointsSize(i);
    k = mod(i-1,6)+1;
    for j = 1:Pointnumber
        testpoint = Points{i,j};
        plot3(testpoint(1), testpoint(2), testpoint(3),symbolchange{k},'MarkerEdgeColor','k',...
            'MarkerFaceColor',growth{i,j},'MarkerSize',10);
    end
end

%Plot lines connecting input points
distance = 0;
for i = 1:Timesteps
    Pointnumber = PointsSize(i);
    k = mod(i-1,6)+1;
    for j = 1:Pointnumber-1
        testpoint = Points{i,j};
        testpoint2 = Points{i,j+1};
        plot3([testpoint(1),testpoint2(1)],[testpoint(2),testpoint2(2)],[testpoint(3),...
            testpoint2(3)],'color',colourchange{k}.*0.7);
    end
end
