function saveBin(varargin)
%
%   saveBin(X,fname)
%   saveBin(X,fname,type)
%   ---------------------
%
%   inputs:
%           X: the data list to be saved
%       fname: pathname+filename
%        type: the type of data to be saved
%              by default, the same as that of input data list
% 	
%   Modified by Yang CHEN 2018.02.24
%            add an option to save multicolomn matrix

X = varargin{1};
fname = varargin{2};

disp(['saving bin-file: ',fname]);

siz = size(X);

if nargin==2
    typeMAT = class(X);    
elseif nargin==3
    typeMAT = varargin{3};
end

% number of points
npts = size(X(:),1);

% data type
typeVTK = typeMAT2VTK(typeMAT);

% write in the file
fid = fopen(fname,'w');

if siz(2)==1
    fprintf(fid, [num2str(npts),'\n']);
    fprintf(fid, [typeVTK,'\n']);
else
    fprintf(fid, [num2str(siz),'\n']);
    fprintf(fid, [typeVTK,'\n']);
end
    
    
fwrite(fid,X(:),typeMAT,0,'b');

%
fclose(fid);

end



function [typeVTK] = typeMAT2VTK(typeMAT)
    
    if (strcmp(typeMAT,'logical')==1)
        typeVTK = 'bit';
    elseif (strcmp(typeMAT,'uint8')==1)
        typeVTK = 'unsigned_char';
    elseif (strcmp(typeMAT,'int8')==1)
        typeVTK = 'char';
    elseif (strcmp(typeMAT,'uint16')==1)
        typeVTK = 'unsigned_short';
    elseif (strcmp(typeMAT,'int16')==1)
        typeVTK = 'short';
    elseif (strcmp(typeMAT,'uint32')==1)
        typeVTK = 'unsigned_int';
    elseif (strcmp(typeMAT,'int32')==1)
        typeVTK = 'int';
    elseif (strcmp(typeMAT,'uint64')==1)
        typeVTK = 'unsigned_long';
    elseif (strcmp(typeMAT,'int64')==1)
        typeVTK = 'long';
    elseif (strcmp(typeMAT,'single')==1)
        typeVTK = 'float';
    elseif (strcmp(typeMAT,'double')==1)
        typeVTK = 'double';
    end
end





