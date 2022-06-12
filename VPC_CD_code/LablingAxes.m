function LablingAxes(SliceDir)
switch SliceDir
    case 'Uz';      xlabel('y (pixel)');    ylabel('x (pixel)');
    case 'Uy';      xlabel('x (pixel)');    ylabel('z (pixel)');
    case 'Ux';      xlabel('z (pixel)');    ylabel('y (pixel)');
end