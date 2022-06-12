function WhitLine(C)
for i=1:size(C,1)
    j=1;
    while sum(C(i,j:j+25)) ~= 0
        j=j+1;
    end
    SavePI(i) = i;
    SavePj(i) = j;
end

hold on; plot(SavePj, SavePI,'-w','LineWidth',3)
            