% find the lentgh of the crack = count(1,:)  
% the y coordiante of the crack = count(2,:)  
% the x coordinate  of the crack = count(3,:)  
function [counter]=findit(cPC_bin)
[Y,X]=size(cPC_bin);
x=1;
counter=zeros(Y,3);

for y=1:Y
    if cPC_bin(y,x)==1
        count=1;
        yy=y;
        for XX=1:X
            if cPC_bin(yy,XX)==1 
                count=count+1;
            else %first corner
                for trial=1:3
                    if cPC_bin(yy+trial,XX)==1
                        yy=yy+trial;
                    elseif cPC_bin(yy-trial,XX)==1
                        yy=yy-trial;
                    end
                end
            end
        end
        counter(y,1)=count; %crack lentgh
        counter(y,2)=yy;    %crack start at y
        counter(y,3)=XX;    
    end
end
