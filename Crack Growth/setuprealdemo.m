%Give pairs of x-y-z coordinates for points at each timestep
%Points(Timestep, point number) = (x,y)
Xdata = {x00000; x00500};
Ydata = {y00000; y00500};
Zdata = {z00000; z00500};

Timepoints = [00000 00500];
Timesteps = size(Timepoints,2);
Points = cell(Timesteps,1);

for i = [1:Timesteps]
    PointsSize(i) = size(Xdata{i},1);
    for j = [1:PointsSize(i)]
        Points{i,j} = [Xdata{i}(j), Ydata{i}(j), Zdata{i}(j)];
    end
end


%De-duplicate
dup = 1; 
while dup == 1
    dup = 0;
    for i = 1:Timesteps-1
        for ii = i:i+1
            PointsSize(i) = size(Xdata{i},1);
            PointsSize(ii) = size(Xdata{ii},1);
            for j = 1:PointsSize(i)
                for k = 1:PointsSize(ii)
                    if isequal(Points{i,j}, Points{ii,k}) & ((j~=k) || (i~=ii))
                        disp (sprintf('Duplicate found. testtime1 %d testtime2 %d j %d k %d', i, ii, j, k ))
                        Points{i,j} = Points{i,j} + [rand rand rand]/100;
                        dup = 1;
                    end
                end
            end
        end
    end
end
    