function frames = readframes(f, GOP_size, start_sequence, DH1)
% Reads Y CIF frames from y4m files
% Read GOP_Size frames from the file with file handle f from 
% frame start_sequence.
% DH1 used to align frames in some sequences e.g. Mobile

H1      = 44;           % y4m header length
H2      = 6;            % y4m header length
CIF     = 352*288;      % Gray CIF frame size in pixels

H1 = H1 + DH1;

frames = cell(1,GOP_size);
index = 1;
for i = start_sequence: start_sequence+GOP_size-1
    CurPos = (i-1)*(H2+CIF+CIF/2)+H1+H2;
    fseek(f,CurPos,-1);
    Img = reshape(fread(f,CIF),352,288)';
    frames{index} = double(Img);
    index = index +1;
end
end

function frames = readframes_sif(f, GOP_size, start_sequence, DH1)

H1      = 44;           % y4m header length
H2      = 6;            % y4m header length
CIF     = 352*240;      % Gray CIF frame size in pixels

H1 = H1 + DH1;

frames = cell(1,GOP_size);
index = 1;
for i = start_sequence: start_sequence+GOP_size-1
    CurPos = (i-1)*(H2+CIF+CIF/2)+H1+H2;
    fseek(f,CurPos,-1);
    Img = reshape(fread(f,CIF),352,240)';
    frames{index} = double(Img);
    index = index +1;
end

end

function PlotRes(RES,IMGS,Set)

figure(2);
for i = 1:size(RES,1)
    
    subplot(2,3,i);
    plot(RES(i,:,1),'.-b');
    hold on;
    plot(RES(i,:,4),'.-r');
    grid off; grid on; grid minor;
    ylabel('PSNR dB'); xlabel('Frame number');
    
    Seq = IMGS{Set(i)}; Seqp = strfind(Seq,'_'); Seq(Seqp) = '-';
    title(Seq);
    
end
sgtitle('PSNR');

figure(3);
for i = 1:size(RES,1)
    
    subplot(2,3,i);
    plot(RES(i,:,2),'.-b');
    hold on;
    plot(RES(i,:,5),'.-r');
    grid off; grid on; grid minor;
    ylabel('MS-SSIM'); xlabel('Frame number');
    
    Seq = IMGS{Set(i)}; Seqp = strfind(Seq,'_'); Seq(Seqp) = '-';
    title(Seq);
    
end
sgtitle('MS-SSIM');

end

function Rx = AL_DCT_ZZ(Img,CR,B)

[H,W]  = size(Img);
[~,zz] = Zz(B);

m  = round(B*B*CR);
Rx = zeros(H,W);

for r = 1:B:H
    for c = 1:B:W
        Patch = Img(r:r+B-1,c:c+B-1);
        D = dct2(Patch);
        D(zz(m+1:end)) = 0;
        Rx(r:r+B-1,c:c+B-1) = idct2(D);
    end
end
end

function Rx = Filter(Ref,Cur,B)

[H,W] = size(Ref);

for r = 1:B:H
    for c = 1:B:W
        Patch = Cur(r:r+B-1,c:c+B-1);
        D = dct2(Patch);
        F = (D~=0);
        Patch = Ref(r:r+B-1,c:c+B-1);
        D = dct2(Patch);
        D = D.*F;
        Rx(r:r+B-1,c:c+B-1) = idct2(D);
    end
end

end

