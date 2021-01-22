function nBTCs = jz_calc_hf_THI(Img,B,CF)
% INPUT
%   Img   = Input Image of size (H,W)
%           H and W should be excatly divisible by B
%   B     = Block size
%   CF    = Compression Factor N/M
%
% OUTPUT
%   nBTCs = Contains the total number of coefficients M1+M2 to collect from each block
%           stored at positions nBTCs(1:B:H,1:B:W)

% Uses the number of the largest M/2 coefficients in a block from Phase1
% to calculate the Phase2 blocks.
% We limit measurements to B^2-1 and fill up blocks to exactly capture M
% measurements.

[H,W] = size(Img);

[~,zz] = Zz(B);

nBTCs = zeros(H,W);

m = round(B*B/CF);
M = round(H*W/CF);


% --------------------------------------------------------
% Phase 1 L-DCT-ZZ with m/2 measurements
TCs = zeros(H,W);
for r = 1:B:H
    for c = 1:B:W
        Patch = Img(r:r+B-1,c:c+B-1);
        D = dct2(Patch);
        TCs(r:r+B-1,c:c+B-1) = D;
    end
end
p = zz(fix(m/2)+1:end);
Mask = ones(B,B);
Mask(p) = 0;
for r = 1:B:H
    for c = 1:B:W
        TCs(r:r+B-1,c:c+B-1) = TCs(r:r+B-1,c:c+B-1).*Mask;
    end
end


% --------------------------------------------------------
% Use Phase 1 measurements to estimate the number of TCs
% to measure in Phase 2.

TCs_sorted = sort(abs(TCs(:)),'descend');
T = TCs_sorted(fix(M/4)+1);

F = abs(TCs)>T;

for r = 1:B:H
    for c = 1:B:W
        nBTCs(r,c) = m/2 + 2*sum(sum(F(r:r+B-1,c:c+B-1)));
        if nBTCs(r,c)>B*B
            nBTCs(r,c) = B*B-1;
        end
    end
end

TCs = floor(nBTCs(1:B:H,1:B:W));
TCs_tot = sum(TCs(:));

% We may not have collected all measurements.
% Collect them all.
p = 0; NBs = numel(TCs);
while M > TCs_tot
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























