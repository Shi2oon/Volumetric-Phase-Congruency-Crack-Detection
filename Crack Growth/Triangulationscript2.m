while finished == 0

    %Test for the shortest next line going forward and backward
    if j < Endj
        lj = norm(Currentline{j+1} - Nextline{k});
        lm = norm(Currentline{Endj-1} - Nextline{Endk});
    else
        lj = 1e6;
        lm = 1e6;
    end

    if k < Endk
        lk = norm(Currentline{j} - Nextline{k+1});
        ln = norm(Currentline{Endj} - Nextline{Endk-1});
    else
        lk = 1e6;
        ln = 1e6;
    end

    
    if lj == 0
        disp (sprintf('Currentline{j} %d %d %d  Currentline{j+1} %d %d %d  Nextline{k} %d %d %d  Nextline{k+1} %d %d %d ', Currentline{j}, Currentline{j+1}, Nextline{k}, Nextline{k+1}))
    end

    %if (i == 1) & (Element == 17)
    %    disp (sprintf('Currentline{j} %d Currentline{j+1} %d Nextline{k} %d Nextline{k+1} %d ', Currentline{j}, Currentline{j+1}, Nextline{k}, Nextline{k+1}))
    %end

    if lj <= min([lk lm ln])
        Element = Element + 1;
        %Connect Currentline{j+1} to Nextline{k}
        Triangles{i, Element} = [0 j j+1 k lj];
        %if i==16
        %    disp (sprintf('Triangles{i, Element} i %d Element %d 0 j %d j+1 %d k %d lj %d', i, Element, j, j+1, k, lj))
        %end
        
        PointB = Currentline{j+1};
        PointX = Nextline{k};

%        if (lj > longest)
%            longest = lj;
%        end

        if lj > norm(growthvector{i,j+1})
            growthvector{i,j+1} = Nextline{k}- Currentline{j+1};
        end
        j = j + 1;

    elseif lk <= min([lj lm ln])
        
        Element = Element + 1;
        %Connect Currentline{j} to Nextline{k+1}
        Triangles{i, Element} = [1 j k+1 k lk];
        %if i==16
        %    disp (sprintf('Triangles{i, Element} i %d Element %d 1 j %d k+1 %d k %d lk %d', i, Element, j, k+1, k, lk))
        %end        
        PointB = Currentline{j};
        PointX = Nextline{k+1};

%        if (lk > longest) 
%            longest = lk;
%        end

        if lk > norm(growthvector{i,j})
            growthvector{i,j} = Nextline{k+1} - Currentline{j};
        end
        k = k + 1;
        
    elseif lm <= min([lj lk ln])
        
        Element = Element + 1;
        %Connect Currentline{Endj-1} to Nextline{Endk}
        Triangles{i, Element} = [0 Endj-1 Endj Endk lm];
        %if i==16
        %    disp (sprintf('Triangles{i, Element} i %d Element %d 0 Endj-1 %d Endj %d Endk %d lm %d', i, Element, Endj-1, Endj, Endk, lm))
        %end
        PointB = Currentline{Endj-1};
        PointX = Nextline{Endk};
        
%        if (lm > longest) 
%            longest = lm;
%        end

        if lm > norm(growthvector{i,Endj-1})
            growthvector{i,Endj-1} = Nextline{Endk} - Currentline{Endj-1};
        end
        Endj = Endj - 1;

    else
        Element = Element + 1;
        %Connect Currentline{Endj} - Nextline{Endk-1}
        Triangles{i, Element} = [1 Endj Endk Endk-1 ln];
        %if i==16
        %    disp (sprintf('Triangles{i, Element} i %d Element %d 1 Endj %d Endk %d Endk-1 %d ln %d', i, Element, Endj, Endk, Endk-1, ln))
        %end
        
        PointB = Currentline{Endj};
        PointX = Nextline{Endk-1};

%        if (ln > longest)
%            longest = ln;
%        end

        if ln > norm(growthvector{i,Endj})
            growthvector{i,Endj} = Nextline{Endk-1} - Currentline{Endj};
        end
        Endk = Endk - 1;

    end

    %Plot edges of triangles
    %plot3([PointB(1),PointX(1)],[PointB(2),PointX(2)],[PointB(3),PointX(3)],'color',colourchange{modi});

    if (j == Endj) && (k == Endk)
        finished = 1;
    end


end