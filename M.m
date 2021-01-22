function y = M(Img,nBTCs,B,zz)
% Computes Adaptive compressive measurements function for IDA
% INPUT
%   Img     of size H,W
%   nBTCs   number of Block TCs in nBTCs(1:H:W,1:B:W)
%   B       block size
%   zz      zigzag scan
% OUTPUT
%   y       measurements

    [H,W] = size(nBTCs);
    Img   = reshape(Img,[H,W]);

    y = [];

    for r = 1:B:H
        for c = 1:B:W
            Patch = Img(r:r+B-1,c:c+B-1);
            D = dct2(Patch);
            y = [y; D(zz(1:nBTCs(r,c)))];
        end
    end
                
end