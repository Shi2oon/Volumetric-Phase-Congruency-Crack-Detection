function [indexMD, indexSD]=Index_para(MainDispl,SliceDir)

switch MainDispl
    case 'Ux'
        indexMD = 1;
    case 'Uy'
        indexMD = 2;
    case 'Uz'
        indexMD = 3;
end

switch SliceDir
    case 'Ux'
        indexSD = 1;
    case 'Uy'
        indexSD = 2;
    case 'Uz'
        indexSD = 3;
end
end
