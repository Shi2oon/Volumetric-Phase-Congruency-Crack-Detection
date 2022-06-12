function chkVair(Vair,VxLst,VOL0,VOLres,method)
%% ========================================================================
%%      check whether the chosen gray value of air is appropriate
%% ========================================================================
%% if it is, 
%%          - the histogram should mostly locates between 0 and 1
%% if it isnot
%%          - great part of the histogram will go outside the interval [0,1]
%%  
%% 15/11/2016
%% Yang Chen
%

over = 0;
while over==0

    Vgap = VOL0 - Vair;

    Vres = VOLres(VxLst);

    Qair = Qair_comput(Vres,Vgap(VxLst),intmax(class(VOLres)),method);

    figure;hist(Qair(:),1000);
    title(['Vair=',num2str(Vair),'; method:',method])
    
    runC = input(['Is the air gray value (',num2str(Vair),') ',...
                                            'satisfactory ? [y]/n '],'s');
    
    if isempty(runC);  runC = 'y';  end
    
    if strcmp(runC,'y')
        over=1;
    else
        Vair = input('User input air gray value: ');
    end
    
end
    

