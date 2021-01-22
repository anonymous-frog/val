function  Img = Mt(y,B,nBTCs,zz)
    
    [H,W] = size(nBTCs);
    Img   = zeros(H,W);
    
    p = 1;
    for r = 1:B:H
        for c = 1:B:W
            dp = nBTCs(r,c);
            D = zeros(B,B);
            D(zz(1:dp)) = y(p:p+dp-1);
            Img(r:r+B-1,c:c+B-1) = idct2(D);
            p = p + dp;
        end
    end
    
    Img = reshape(Img,[H*W,1]);
    
end