function nBTCs = jz_calc_hf_3(Img,B,CF)
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
%
% Uses the number of the largest M coefficients in a block from all TCs
% to calculate the Phase2 block measurements.
%
% THIS IS FEASIBLE IFF A PREVIOUS FRAME IS USED EITHER Key or nonKey
%
% We limit measurements to B^2-1 and fill up blocks to exactly capture M
% measurements.


[H,W] = size(Img);

[~,zz] = Zz(B);

nBTCs = zeros(H,W);

m = fix(B*B/CF);
M = fix(H*W/CF);

TCs = zeros(H,W);

for r = 1:B:H
    for c = 1:B:W
        Patch = Img(r:r+B-1,c:c+B-1);
        D = dct2(Patch);
        TCs(r:r+B-1,c:c+B-1) = D;
    end
end

TCs_sorted = sort(abs(TCs(:)),'descend');
T = TCs_sorted(M+1);

F = abs(TCs)>T;

for r = 1:B:H
    for c = 1:B:W
        nBTCs(r,c) = sum(sum(F(r:r+B-1,c:c+B-1)));
    end
end























