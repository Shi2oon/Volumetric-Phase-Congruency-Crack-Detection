function [X] = lect_raw(varargin)
% LECTURE DE FICHIER  *.raw
% =================================
%
%   [X] = lect_raw(nfich0,nx,ny,nz)
%   [X] = lect_raw(nfich0,nx,ny,nz,type)
%   [X] = lect_raw(nfich0,nx,ny,nz,nblocs)
%   [X] = lect_raw(nfich0,nx,ny,nz,type,endianness)
%   [X] = lect_raw(nfich0,nx,ny,nz,type,nblocs)
%   [X] = lect_raw(nfich0,nx,ny,nz,type,endianness,nblocs)
%   [X] = lect_raw(nfich0,nx,ny,nz,type,endianness,nblocs,[n1,n2])
%
% Inputs: 
%   - nfich0  : raw file name
%   - nx,ny,nz: image dimension
%   - type    : (by default 'uint8') the datatype of the volume data to be
%               read 
%   - endianness: machine format
%                 'n' - Your system byte ordering
%                 'b' - Big-endian ordering
%                 'l' - Little-endian ordering
%                 's' - Big-endian ordering, 64-bit long data type
%                 'a' - Little-endian ordering, 64-bit long data type
%   - nblocs  : (by default 1) it is possible to read a huge voluem by
%                nblocs times 
%   - [n1,n2] : the starting and the ending slice numbers
%
% Outputs:
%   - X: the reading result
% ==================================================
%
% Yang CHEN version 2018.12.18

nfich0 = varargin{1};
nx = varargin{2};
ny = varargin{3};
nz = varargin{4};
type = 'uint8';
endianness = 'b';
nblocs = 1;
switch nargin
    case 5
        if ischar(varargin{5})
            type=varargin{5};
        elseif isnumeric(varargin{5})
            nblocs = varargin{5};
        end
    case 6
        type=varargin{5};
        if ischar(varargin{6})
            endianness=varargin{6};
        elseif isnumeric(varargin{6})
            nblocs = varargin{6};
        end
    case 7
        type=varargin{5};
        endianness=varargin{6};
        nblocs=varargin{7};
    case 8
        type=varargin{5};
        endianness=varargin{6};
        nblocs=varargin{7};
        n1 = varargin{8}(1);
        n2 = varargin{8}(2);
end

tic
disp(['reading the raw image file: ',nfich0]);
fid=fopen(nfich0,'r');

if nargin<8 %if the starting/ending slice numbers are not specified
    X = zeros(nx,ny,nz,type);
    inc=ceil(nz/nblocs);
    over = 0;
    for ibloc=1:nblocs
        i1 = 1+inc*(ibloc-1);
        i2 = i1+inc-1;
        if i2>nz; i2=nz; inc=i2-i1+1; over=1; end
        B=fread(fid,nx*ny*inc,type,0,endianness);
        X(:,:,i1:i2)=permute(reshape(B,ny,nx,inc),[2 1 3]);%"permute" to adjust the fread with the conventional coordinates (x<=>column, y<=>row)
        disp(num2str([ibloc,i1,i2,inc]))
        if over==1; break; end
    end

elseif nargin==8 %if the starting/ending slice numbers are specified
    nz = n2-n1+1;
    % move to the the begining slice (n1-th)
    nbytes = (n1-1) * (nx*ny) * nBytes(type);

    fseek(fid,nbytes,'cof');
    
    X = zeros(nx,ny,nz,type);
    inc=ceil(nz/nblocs);  over=0;
    for ibloc=1:nblocs
        i1 = 1+inc*(ibloc-1);
        i2 = i1+inc-1;
        if i2>nz; i2=nz; inc=i2-i1+1; over=1;  end
        B=fread(fid,nx*ny*inc,type,0,endianness);
        X(:,:,i1:i2)=permute(reshape(B,ny,nx,inc),[2 1 3]);%"permute" to adjust the fread with the conventional coordinates (x<=>column, y<=>row)
        disp(num2str([i1+n1-1,i2+n1-1,inc]))
        if over==1; break; end
    end

end
fclose(fid);

toc


function [nbytes] = nBytes(typeMAT)
if (strcmp(typeMAT,'uint8')==1)
    nbytes = 1;
elseif (strcmp(typeMAT,'int8')==1)
    nbytes = 1;
elseif (strcmp(typeMAT,'uint16')==1)
    nbytes = 2;
elseif (strcmp(typeMAT,'int16')==1)
    nbytes = 2;
elseif (strcmp(typeMAT,'uint32')==1)
    nbytes = 4;
elseif (strcmp(typeMAT,'int32')==1)
    nbytes = 4;
elseif (strcmp(typeMAT,'uint64')==1)
    nbytes = 8;
elseif (strcmp(typeMAT,'int64')==1)
    nbytes = 8;
elseif (strcmp(typeMAT,'single')==1)
    nbytes = 4;
elseif (strcmp(typeMAT,'float')==1)
    nbytes = 4;
elseif (strcmp(typeMAT,'double')==1)
    nbytes = 8;
end
end


end