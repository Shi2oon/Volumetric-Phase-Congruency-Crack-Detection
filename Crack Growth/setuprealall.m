%Give pairs of x-y-z coordinates for points at each timestep
%Points(Timestep, point number) = (x,y)
function [Points,PointsSize]=setuprealall(Xdata,Ydata,Zdata,Timepoints)
% Xdata = {x00000; x00500; x01000; x02000; x02500; x03000; x03500; x04000; x05000; x05250; x05500; x05750; x06000; x06250; x06500; x06750; x07000; x07250; x07500; x07750; x08250; x08750; x09250; x10500};
% Ydata = {y00000; y00500; y01000; y02000; y02500; y03000; y03500; y04000; y05000; y05250; y05500; y05750; y06000; y06250; y06500; y06750; y07000; y07250; y07500; y07750; y08250; y08750; y09250; y10500};
% Zdata = {z00000; z00500; z01000; z02000; z02500; z03000; z03500; z04000; z05000; z05250; z05500; z05750; z06000; z06250; z06500; z06750; z07000; z07250; z07500; z07750; z08250; z08750; z09250; z10500};
% Timepoints = [00000 00500 01000 02000 02500 03000 03500 04000 05000 05250 05500 05750 06000 06250 06500 06750 07000 07250 07500 07750 08250 08750 09250 10500];
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
                        Points{i,j} = Points{i,j} + [rand rand rand]/5000;
                        dup = 1;
                    end
                end
            end
        end
    end
end
    
