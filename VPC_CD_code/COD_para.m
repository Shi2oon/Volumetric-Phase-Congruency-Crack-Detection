function [cod,coda] = COD_para(STEP,lastindex,indexMD);
[sizem, sizen] = size(STEP);
if lastindex == 1
    for i=1:sizem
        if isempty(find(isnan(STEP(i,:))==1, 1))==0
            PLEB = find(isnan(STEP(i,:))==1);
            UPPEROP(i) = max(PLEB);
            LOWEROP(i) = min(PLEB);
            X1n = LOWEROP(i)-1;
            X2n = UPPEROP(i)+1;
            X1na = X1n-1;
            X2na = X2n+1;
            X1nb = X1na-1;
            X2nb = X2na+1;
            codW1(i) = abs(STEP(i,X1n)-STEP(i,X2n));
            codW2(i) = abs(mean([STEP(i,X1n), STEP(i,X1na), STEP(i,X1nb)])-mean([STEP(i,X2n), STEP(i,X2na), STEP(i,X2nb)]));
        end
    end
elseif lastindex == 2 && indexMD == 1
    for i=1:sizen
        if isempty(find(isnan(STEP(:,i))==1, 1))==0
            PLEB = find(isnan(STEP(:,i))==1);
            UPPEROP(i) = max(PLEB);
            LOWEROP(i) = min(PLEB);
            X1n = LOWEROP(i)-1;
            X2n = UPPEROP(i)+1;
            X1na = X1n-1;
            X2na = X2n+1;
            X1nb = X1na-1;
            X2nb = X2na+1;
            codW1(i) = abs(STEP(X1n,i)-STEP(X2n,i));
            codW2(i) = abs(mean([STEP(X1n,i), STEP(X1na,i), STEP(X1nb,i)])-mean([STEP(X2n,i), STEP(X2na,i), STEP(X2nb,i)]));
        end
    end
elseif lastindex == 3 && indexMD == 2
    for i=1:sizen
        if isempty(find(isnan(STEP(:,i))==1, 1))==0
            PLEB = find(isnan(STEP(:,i))==1);
            UPPEROP(i) = max(PLEB);
            LOWEROP(i) = min(PLEB);
            X1n = LOWEROP(i)-1;
            X2n = UPPEROP(i)+1;
            X1na = X1n-1;
            X2na = X2n+1;
            X1nb = X1na-1;
            X2nb = X2na+1;
            codW1(i) = abs(STEP(X1n,i)-STEP(X2n,i));
            codW2(i) = abs(mean([STEP(X1n,i), STEP(X1na,i), STEP(X1nb,i)])-mean([STEP(X2n,i), STEP(X2na,i), STEP(X2nb,i)]));
        end
    end
else
    for i=1:sizem
        if isempty(find(isnan(STEP(i,:))==1, 1))==0
            PLEB = find(isnan(STEP(i,:))==1);
            UPPEROP(i) = max(PLEB);
            LOWEROP(i) = min(PLEB);
            X1n = LOWEROP(i)-1;
            X2n = UPPEROP(i)+1;
            X1na = X1n-1;
            X2na = X2n+1;
            X1nb = X1na-1;
            X2nb = X2na+1;
            codW1(i) = abs(STEP(i,X1n)-STEP(i,X2n));
            codW2(i) = abs(mean([STEP(i,X1n), STEP(i,X1na), STEP(i,X1nb)])-mean([STEP(i,X2n), STEP(i,X2na), STEP(i,X2nb)]));
        end
    end
end
cod = codW1;
coda= codW2;
end