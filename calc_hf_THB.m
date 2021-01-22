function [nBTCs] = jz_calc_hf_1(Img,B,CF)
% INPUT
%   Img   = Input Image of size (H,W)
%           H and W should be excatly divisible by B
%   B     = Block size
%   CF    = Compression Factor N/M
%   ADAPT = Adaptation algorithm. ADAPT = 0 > no adaptation
%
% OUTPUT
%   nBTCs = Contains the total number of coefficients M1+M2 to collect from each block
%           stored at positions nBTCs(1:B:H,1:B:W)


[H,W] = size(Img);

[~,zz] = Zz(B);

DH = zeros(H,W);
nBTCs = zeros(H,W);     % Will contain numbger of TCs=M_1+M_2 in nBTCs(1:B:H,1:B:W)
nTCs = fix(H*W/CF);     % Number of transfer coefficients = M
nB = (H*W)/(B*B);       % Number of blocks
nTCsB = fix(nTCs/nB);   % Number TCs per block = Mb. nTCsB/2 = M1

% T optimized empirically on Set_final_{256|512} - previously T = 15;

if H >= 512
    if CF >= 10
        T = 60;
    elseif CF <= 2
        T = 15;
    else
        T = 30;
    end
elseif H <= 256
    if CF >= 1/0.3
        T = 60;
    elseif CF <= 2
        T = 15;
    else
        T = 30;
    end
else
    T = 30;
end

F = 0.50; % Fraction of TCs collected initially

for r = 1:B:H
    for c = 1:B:W
        Patch = Img(r:r+B-1,c:c+B-1);
        D = dct2(Patch);
        Dz = D(zz);
        Dz(fix(nTCsB+1):end)=0;
        DH(r,c) = sum(abs(Dz(1:int16(nTCsB*F)))>T);  % 33.57 dB
    end
end

DHT = sum(abs(DH(:)));

for r = 1:B:H
    for c = 1:B:W
                        
        % nBTCs(r,c) = M2 + M1 and nTCsB/2 = M1
        nBTCs(r,c) = round((nTCs-nTCsB*F*nB)*DH(r,c)/DHT)+nTCsB*F;
        
        if nBTCs(r,c)>(B^2-1)
            nBTCs(r,c)=B^2-1;
        end
                
    end
end

TCs = floor(nBTCs(1:B:H,1:B:W));
TCs_tot = sum(TCs(:));

%We may not have collected enough measurements
p = 0; NBs = numel(TCs);
while nTCs > TCs_tot
    p = p + 1;
    if p > NBs
        p = 1;
    end
    if TCs(p)<(B*B-1)
        TCs(p) = TCs(p)+1;
        TCs_tot = sum(TCs(:));
    end
end

TCs = uint32(TCs);
nBTCs = zeros(H,W);
nBTCs(1:B:H,1:B:W) = TCs;

end






















