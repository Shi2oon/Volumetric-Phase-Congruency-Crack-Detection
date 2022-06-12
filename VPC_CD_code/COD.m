function [cod,coda,X2c,Y2c,KY,KX,CPX,CPY] = COD(STEP);

[sizem, sizen] = size(STEP);
    for i=1:sizem
        if isempty(find(isnan(STEP(i,:))==1, 1))==0
            PLEB = find(isnan(STEP(i,:))==1);
            UPPEROP(i) = max(PLEB);
            LOWEROP(i) = min(PLEB);
    %         X1 = 1:LOWEROP(i)-1;
    %         X2 = UPPEROP(i)+1:sizen

            X1n = LOWEROP(i)-1;
            X2n = UPPEROP(i)+1;
            X1na = X1n-1;
            X2na = X2n+1;
            X1nb = X1na-1;
            X2nb = X2na+1;
            codW1(i) = abs(STEP(i,X1n)-STEP(i,X2n));
            codW2(i) = abs(mean([STEP(i,X1n), STEP(i,X1na), STEP(i,X1nb)])-mean([STEP(i,X2n), STEP(i,X2na), STEP(i,X2nb)]));
%             codW1(i) = abs(STEP(X1n,i)-STEP(X2n,i));
%             codW2(i) = abs(mean([STEP(X1n,i), STEP(X1na,i), STEP(X1nb,i)])-mean([STEP(X2n,i), STEP(X2na,i), STEP(X2nb,i)]));

    %         codW2(i) = abs(mean([STEP(i,X1na), STEP(i,X1nb)])-mean([STEP(i,X2na), STEP(i,X2nb)]));
%             X1a(i) = X1STEP(i,X1n);
%             X1b(i) = X1STEP(i,X2n);
%             Y1a(i) = Y1STEP(i,X1n);
%             Y1b(i) = Y1STEP(i,X2n);

     
        
        
        
        %         pU = polyfit(X1,STEP(i,X1),1);
%         pD = polyfit(X2,STEP(i,X2),1);
%         pUA = mean(STEP(i,X1));
%         pDA = mean(STEP(i,X2));
%         p1(i,:) = [pU(1) pD(1)];
%         p2(i,:) = [pU(2) pD(2)];
%         cod1(i) = abs(pU(2) - pD(2));
%         codT(i) = abs(pUA - pDA);
%         MODEL(i,X1) = X1*pU(1) + pU(2);
%         MODEL(i,X2) = X2*pD(1) + pD(2);
%         MODEL1(i,X1) = pUA;
%         MODEL1(i,X2) = pDA;
%         if STEP(i,(X2(end)))>STEP(i,(X1(1)))
%             codW1(i) = STEP(i,(X2(end))) - STEP(i,(X1(1)));
%         else codW1(i) = STEP(i,(X1(end))) - STEP(i,(X2(1)));
%             
%         end



%     else
%         X3 = 1:sizen;
%         p = polyfit(X3,STEP(i,X3),1);
%         pA = mean(STEP(i,X3));
%         MODEL(i,X3) = X3*p(1)+p(2);
%         MODEL1(i,X3) = pA;
%     end
    end
    
    

    
end

% 
% 
% 
% % cod = cod1;
% % codTA = codT;
% % codTA(end+1) = 0;
% % cod(end+1) = 0;
% % FMODEL = MODEL;
% % MODEL1(MODEL1==0)=nan;
% % MMODEL = MODEL1;
% % codW1(end+1) = 0;
% %     X1a(X1a==0) = [];
% %     X1b(X1b==0) = [];
% %     Y1a(Y1a==0) = [];
% %     Y1b(Y1b==0) = [];
% %     length(X1a)
% %     length(X1b)
% %     length(Y1a)
% %     length(Y1b)
% CPX = (X1a + X1b)/2;
% figure
% plot(CPX)
% hold on
% plot(X1a,'r')
% plot(X1b,'c')
% plot(1:length(Y1a),1024,'-b')
% KY = Y1a;
% CPY = (Y1a + Y1b)/2;
% figure
% plot(CPY)
% hold on
% plot(1:length(Y1a),1024,'-b')
% plot(Y1a,'r')
% plot(Y1b,'c')
% KX = X1a;
%     X2c =[X1a, X1b];
%     Y2c =[Y1a, Y1b];
%     X2c(X2c==0) = [];
%     Y2c(Y2c==0) = [];
    cod = codW1;
    coda= codW2;
%     cod(cod==0) = [];
%     coda(coda==0) = [];
%     CPX(CPX==0) = [];
%     CPY(CPY==0) = [];