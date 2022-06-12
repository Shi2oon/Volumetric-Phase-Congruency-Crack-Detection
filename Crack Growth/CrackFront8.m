function Data = CrackFront8(Points, PointsSize, Timepoints, SamplingSteps)

Timesteps  = size(Points,1);

%Load colormap
load 'mycolormap.mat'

disp (sprintf('%d timesteps detected \n', Timesteps ))

handlefirstplot = figure;
hold on
daspect([1 1 1])
title('Nodes inserted into intial mesh')

%lighting gouraud
%material dull   
%light('Position',[1 0 0],'Style','infinite');
%light('Position',[0 1 0],'Style','infinite');
%light('Position',[0 0 1],'Style','infinite');

colourchange = {[1,0,0] [0.5,0.5,0] [0,1,0] [0,1,1] [0,0,1] [1,0,1]};
symbolchange = {'o' 's' '^' '<' 's' 'o'};

%Plot lines connecting input points before resampling, and calculate average point spacing

FirstPlots

disp('PointsSize before resampling ');
disp(sprintf( ' %d', PointsSize));

%Resample long lines to add extra points and plot smaller points for
%inserted points
if SamplingSteps > 0 
    ResampleScript

    Points = InterPoints;
    PointsSize = InterPointsSize;
    disp('PointsSize after resampling ');
    disp(sprintf( ' %d', PointsSize));
end

%De-duplicate again in case we've created duplicates during resampling
dup = 1; 
while dup == 1
    dup = 0;
    for i = 1:Timesteps-1
        for ii = i:i+1
            for j = 1:PointsSize(i)
                for k = 1:PointsSize(ii)
                    if isequal(Points{i,j}, Points{ii,k}) && ((j~=k) || (i~=ii))
                        disp (sprintf('Duplicate found after resampling. testtime1 %d testtime2 %d j %d k %d', i, ii, j, k ))
                        Points{i,j} = Points{i,j} + [rand rand rand]/5000;
                        dup = 1;
                    end
                end
            end
        end
    end
end

MaxPointnumber = max(PointsSize);
growthvector = cell(Timesteps-1,MaxPointnumber);

%Define triangles connecting points.

%Loop to find nearest point on next line
NearestPointC = cell(Timesteps-1,MaxPointnumber);
for i = 1:Timesteps-1
    
    % Set up Currentline and NextLine Arrays
    Pointnumber = PointsSize(i);
    PointnumberNextrow = PointsSize(i+1);

    Currentline = cell(Pointnumber);
    for j = 1:Pointnumber
        Currentline{j} = Points{i,j};
    end

    Nextline = cell(PointnumberNextrow);
    for k = 1:PointnumberNextrow
        Nextline{k} = Points{i+1,k};
    end

    %Find nearest point on next line
    for j = 1:Pointnumber
        mindist = 1e6;
        for k = 1:PointnumberNextrow
            dist = norm(Currentline{j} - Nextline{k});
            if dist < mindist
                mindist = dist;
                NearestPointC{i,j}(1) = k;
            end
        end
    end
end


%Loop to find nearest points on previous line
NearestPointPrevious = cell(Timesteps-1,MaxPointnumber);
for i = 1:Timesteps-1
    
    % Set up Currentline and NextLine Arrays
    Pointnumber = PointsSize(i);
    PointnumberNextrow = PointsSize(i+1);

%    disp (sprintf('Timestep %d Pointnumber %d PointnumberNextrow %d \n', i, Pointnumber, PointnumberNextrow  ))

    modi = mod(i-1,6)+1;
    
    for j = 1:Pointnumber
        Currentline{j} = Points{i,j};
    end

    for k = 1:PointnumberNextrow
        Nextline{k} = Points{i+1,k};
    end

    %Find nearest point on previous line
    for k = 1:PointnumberNextrow
        mindist = 1e6;
        for j = 1:Pointnumber
            dist = norm(Currentline{j} - Nextline{k});
            if dist < mindist
                mindist = dist;
                NearestPointPrevious{i+1,k}(1) = j;
            end
        end
    end
end

%for i = 1:Timesteps-1
%    disp (sprintf('i    %d', i ))
%    
%    Pointnumber = PointsSize(i);
%    for j = 1:Pointnumber
%        disp (sprintf('j %d C %d', j, NearestPointC{i,j}(1)  ))
%    end
%    
%    PointnumberNextrow = PointsSize(i+1);
%    for j = 1:PointnumberNextrow
%        disp (sprintf('j %d C %d', j, NearestPointPrevious{i+1,j}(1)  ))
%    end
%end


%Test which points mirror each other in nearest neighbour matrix
%Set up regions to triangulate between lines.

ConnectionsNo = zeros(Timesteps-1);
Connections = zeros(Timesteps-1,MaxPointnumber);

for i = 1:Timesteps-1
    Pointnumber = PointsSize(i);
    PointnumberNextrow = PointsSize(i+1);
    m = 1; %m is the counter for lines to use for triangulations

    klast = 0;
    for j = 1:Pointnumber
        k = NearestPointC{i,j}(1);
        l = NearestPointPrevious{i+1,k}(1);
%        disp (sprintf('Timestep %d j %d k %d l %d', i, j, k, l ))

        if  (j == l) && (k <= klast)
            %This plots larger markers on out of order points.
            %plot3(Points{i,j}(1), Points{i,j}(2), Points{i,j}(3),symbolchange{modi},'MarkerEdgeColor','k',...
            %'MarkerFaceColor',colourchange{modi},'MarkerSize',20); 
            %disp (sprintf('Out of order connection. i %d j %d k %d l %d', i, j, k, l ))
            NearestPointC{i,j}(2) = 0;
        else
            if (j == l)
                 NearestPointC{i,j}(2) = 1;
                 Connections(i,m) = j;
                 %disp (sprintf('Connection found. i %d j %d k %d l %d m %d', i, j, k, l, m ))
                 m = m + 1;
                 klast = k;
            else 
                 NearestPointC{i,j}(2) = 0;
            end
        end
    end
    ConnectionsNo(i) = m - 1;
end


%Loop to perform triangulation

%longest  = 0;
NElements = zeros(Timesteps-1);
for i = 1:Timesteps-1
    Stepsize = Timepoints(i+1) - Timepoints(i);

    if i == 16
        disp (sprintf('Timestep %d Pointnumber %d PointnumberNextrow %d ', i, Pointnumber, PointnumberNextrow  ))
    end
    modi = mod(i-1,6)+1;
    Pointnumber = PointsSize(i);
    PointnumberNextrow = PointsSize(i+1);
    
    for j = 1:Pointnumber
        Currentline{j} = Points{i,j};
    end

    for j = 1:PointnumberNextrow
        Nextline{j} = Points{i+1,j};
    end
    

    growthvector{i,1} = Nextline{1}- Currentline{1};
%    if (norm(growthvector{i,1})/Stepsize > longest)
%        longest = norm(growthvector{i,1})/Stepsize;
%    end

    Element = 0;

    
    for m = 0:ConnectionsNo(i)
        
        if m == 0
            j = 1;
            k = 1;
        else
            j = Connections(i,m);
            k = NearestPointC{i,j}(1);
        end
        
        if m ~= ConnectionsNo(i)
            Endj = Connections(i,m+1);
            Endk = NearestPointC{i,Endj}(1);
        else
            Endj = Pointnumber;
            Endk = PointnumberNextrow;
        end
        
        %if i==16
        %    disp (sprintf('Triangulating i %d j %d k %d  m %d Endj %d Endk %d', i, j, k, m, Endj, Endk))
        %end

        if (j == Endj) && (k == Endk)
            finished = 1;
        else
            finished = 0;
        end

        li = norm(Currentline{j} - Nextline{k});
        if li > norm(growthvector{i,j})
            growthvector{i,j} = Nextline{k}- Currentline{j};
%            if (li/Stepsize > longest) 
%                longest = li/Stepsize;
%            end
        end

        lo = norm(Currentline{Endj} - Nextline{Endk});
        if lo > norm(growthvector{i,Endj})
            growthvector{i,Endj} = Nextline{Endk}- Currentline{Endj};
 %           if (lo/Stepsize > longest) 
 %               longest = lo/Stepsize;
 %           end
        end
        
        Triangulationscript2
    end
    
    NElements(i) = Element;
    %disp (sprintf('NElements(i) %d \n', NElements(i) ))

end

%disp (sprintf('Longest connection %d', longest))

%Plot  figure with vectors drawn in for longest connections.
handlevectors = figure;
hold on
daspect([1 1 1])
title('Triangulation of nodes')

%Call firstPlots to plot all points and calculate new average point spacing
FirstPlots

%Plot connections to next line
for i = 1:Timesteps-1
    Pointnumber = PointsSize(i);
    modi = mod(i-1,6)+1;
    for j = 1:Pointnumber
        if NearestPointC{i,j}(2) == 1
            %Plot nearest point on next line in different line style
            PointB = Points{i,j};
            PointX = Points{i+1, NearestPointC{i,j}(1)};
            plot3([PointB(1),PointX(1)],[PointB(2),PointX(2)],[PointB(3),PointX(3)],'color',[0.5 0.5 0.5],'linewidth', 3);
        else
            PointB = Points{i,j};
            PointX = Points{i,j} + growthvector{i,j};
            plot3([PointB(1),PointX(1)],[PointB(2),PointX(2)],[PointB(3),PointX(3)],'color',colourchange{modi});
        end
    end
end

%Calculate the normalised velocity forward from each point
k = 0;
for i = 1:Timesteps-1
    Pointnumber = PointsSize(i);
    Stepsize = Timepoints(i+1) - Timepoints(i);

    for j = 1:Pointnumber
        k = k+1;

        g = norm(growthvector{i,j});
        growth1d(k) = g/Stepsize;
        growth2d(i,j) = g/Stepsize;
        length2d(i,j) = g;

        g1 = log2(g)/2 + 2; %red
        g2 = log2(g)/2 + 4; %green
        g3 = log2(g)/2 + 5; %blue
        if g1 > 1
            g1 = 3-g1;
        end
        if g1 > 1
            g1 = 1;
        end

        if g1 < 0
            g1 = 0;
        end

        if g2 > 1
            g2 = 4-g2;
        end
        if g2 > 1
            g2 = 1;
        end
        if g2 < 0
            g2 = 0;
        end

        if g3 > 1
            g3 = 3 - g3;
        end
        if g3 > 1
            g3 = 1;
        end
        if g3 < 0
            g3 = 0;
        end

        growth{i,j} = [g1 g2 g3];
    end

    meangrowth(i) = mean(growth2d(i,:));
    meanlength(i) = mean(length2d(i,:));

%Enable this to plot histograms at each point    
%handlehiststep(i) = figure;
%hist(growth2d(i,:),[0:0.01:1])

end

fastestperstep = max(growth2d(:,:)');
fastest = max(fastestperstep);
longest = max(length2d(:,:)');

Totalpoints = k;

for i = 1:Timesteps-1
    Pointnumber = PointsSize(i);
    Stepsize = Timepoints(i+1) - Timepoints(i);

    for j = 1:Pointnumber
        growth2dnorm(i,j) = growth2d(i,j)/fastest;
    end
end


handlemean = figure;
plot(Timepoints(1:Timesteps-1),meangrowth)
title('Mean growth rates')
xlabel('Time')
ylabel('Growth rates')

handlelength = figure;
plot(Timepoints(1:Timesteps-1),meanlength)
title('Mean lengths')
xlabel('Time')
ylabel('Length')

%Plot the triangles with the colourmap set up properly and set
%growth vector

handlesmoothplot = figure;
hold on
daspect([1 1 1])
colormap(mycolormap)
title('Surface showing normalised crack velcoity')

for i = 1:Timesteps-1
    for j = 1 : NElements(i)
        node = Triangles{i, j};

        if numel(node) == 0
            disp (sprintf('No nodes for Triangle. i %d j %d Triangles(i,j-1) %d %d %d %d', i, j, Triangles{i, j-1}(1),Triangles{i, j-1}(2),Triangles{i, j-1}(3),Triangles{i, j-1}(4) ))
        else
            X = [Points{i,node(2)}(1) Points{i+node(1),node(3)}(1) Points{i+1,node(4)}(1)];
            Y = [Points{i,node(2)}(2) Points{i+node(1),node(3)}(2) Points{i+1,node(4)}(2)];
            Z = [Points{i,node(2)}(3) Points{i+node(1),node(3)}(3) Points{i+1,node(4)}(3)];

            if i == Timesteps-1
                patch(X,Y,Z, [ growth2dnorm(i,node(2)); growth2dnorm(i,node(2)); growth2dnorm(i,node(2)) ] )
            else
                patch(X,Y,Z, [ growth2dnorm(i,node(2)); growth2dnorm(i+node(1),node(3)); growth2dnorm(i+1,node(4)) ] )
            end
        end
    end
end

%Plot figure with coloured points to represent the velocity of the crack
%velocitiesfigurescript

handlehist = figure;

[n,xout] = hist(growth1d/fastest,[0:0.004:1]);
bar(xout,n/Totalpoints*100);
title('Histogram of normalised growth rates over all timesteps')
xlabel('Normalised crack growth rate')
ylabel('Count as percentage')

Data = {Points, growthvector, growth2d, length2d, growth2dnorm, longest, fastest, fastestperstep, meangrowth, meanlength};

end