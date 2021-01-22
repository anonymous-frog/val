function [zz,z] = ZZ(B)
% JPEG zigzag scan function
% B is the blocksize

zz = zeros(B,B);
z  = zeros(1,B*B);

zz(1,1)=1;
z(1,1)=1;

for c = 2:2:B
    n = c*(c-1)/2+1;
    r = 1;
    while c~=0
        zz(r,c) = n;
        z(n)=(sub2ind([B,B],r,c));
        c = c - 1;
        r = r + 1;
        n = n + 1;
    end
end

N = n;

for r = 3:2:B
    n = r*(r-1)/2+1;
    c = 1;
    while r~=0
        zz(r,c) = n;
        z(n)=(sub2ind([B,B],r,c));
        c = c + 1;
        r = r - 1;
        n = n + 1;
    end
end

n = N;

for c = 2:2:B
    r = B;
    C = c;
    while C<(B+1)
        zz(r,C) = n;
        z(n)=(sub2ind([B,B],r,C));
        C = C + 1;
        r = r - 1;
        n = n + 1;
    end
    n = n + B - c;
end

n = N + B - 1;

for r = 3:2:B
    c = B;
    R = r;
    while R<(B+1)
        zz(R,c) = n;
        z(n)=(sub2ind([B,B],R,c));
        c = c - 1;
        R = R + 1;
        n = n + 1;
    end
    n = n + B - r;
end

z = z';

end


