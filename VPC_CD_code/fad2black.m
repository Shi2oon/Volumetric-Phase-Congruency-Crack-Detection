function [visPC,crackLength]=fad2black(data,spacing,output)
    im1 = -data;
        
    % fill in NaNs
    im2 = inpaint_nans(im1,1);
    
    % outlier deletion
    [im3] = OutDel(im2,10,0.2);
    im4 = inpaint_nans(im3,1);
    
    % phase congruency
    PC = phasecongmono(im4, 100, 3, 2.5, 0.55, 2.0);
    
    % crop PC
    %cPC = PC(100:150,:);
    cPC = PC;
    % binarize (threshold) image
    cPC_bin = imbinarize(cPC,0.004); 
    %imagesc(cPC_bin)
    %contourf(data.xMap_px,data.yMap_px,cPC_bin)
    
    % find the lentgh of the crack = count(1,:)  
    % the y coordiante of the crack = count(2,:)  
    % the x coordinate  of the crack = count(3,:) 
   
    
    % calculate crack length from size of largest region
    cc = bwconncomp(cPC_bin);
    rp = regionprops(cc);
    rp = struct2cell(rp);
    rp = rp';
    rp = sortrows(rp, 1);
    boundingbox = rp{end,3};
    crackLength = boundingbox(1)+boundingbox(3);
    crackLength = crackLength * spacing;

%     fprintf('Crack length = %3.0f pixels\n',crackLength)
    
    visPC.im1 = im1;
    visPC.im2 = im2;
    visPC.im3 = im3;
    visPC.im4 = im4;
    visPC.cPC = cPC;
    visPC.cPC_bin = cPC_bin;
    
drawPCVis( visPC, iStage,output);