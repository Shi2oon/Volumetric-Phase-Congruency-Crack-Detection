function saveImRAW(fnameS,V)
%save image to raw data
%
%   saveImRAW(fnameS,V)
%   -------------------
%
%   Inputs:
%       - fnameS : file name
%       -      V : image to be saved
%
% Yang CHEN 2019.02.03
%
    

fprintf(['saving the image to ',fnameS,'\n']);

% file information
volinfo = ['width = ',num2str(size(V,1)),'    '...
           'height = ',num2str(size(V,2)), '    '...
           'number of slices = ',num2str(size(V,3)), '    ',...
           class(V)];
fnameSinfo = [fnameS,'info'];
fid = fopen(fnameSinfo,'w');
fprintf(fid,'%s',volinfo);
fclose(fid);

% file saving
fid = fopen(fnameS,'w');
fwrite(fid,V(:),class(V));
fclose(fid);

fprintf('saving complete\n');
