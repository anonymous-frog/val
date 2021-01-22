function [x_hat,nBTCs,Time,CFa,PSNR] = Codec_Type(Img,FrameType,AdaptType,B,ADAPT,CF,CS,iters,DF,denoiser)

    persistent Ref1       % first reference frame
    persistent Ref2       % second reference frame
    persistent PFrame     % previous frame

    [H,W] = size(Img);
    PSNR = [];
    [~,zz] = Zz(B);

    if ADAPT == 0
        nBTCs = zeros(H,W);
        m = fix(B*B/CF);
        for r = 1:B:H
            for c = 1:B:W
                nBTCs(r,c) = m;
            end
        end
        CFa = (H*W)/(sum(nBTCs(:)));
    else
        switch AdaptType
            case 1 % THB-DD
                nBTCs = fix(calc_hf_THB(Img,B,CF));
                CFa = (H*W)/sum(nBTCs(:));
            case 2 % THI-DD
                nBTCs = fix(calc_hf_THI(Img,B,CF));
                CFa = (H*W)/sum(nBTCs(:));
            case 3 % Ref - TDD
                nBTCs = fix(calc_hf_TDD(Ref1,B,CF));
                CFa = (H*W)/sum(nBTCs(:));
            case 4 % PFrame - TDD
                nBTCs = fix(calc_hf_TDD(PFrame,B,CF));
                CFa = (H*W)/sum(nBTCs(:));
            case 5 % PFrame - TDD
                nBTCs = fix(calc_hf_MTDD(Img,Ref1,B,CF));
                CFa = (H*W)/sum(nBTCs(:));
        end
    end
    
    y = M(Img,nBTCs,B,zz);
    tic
    x_hat = Mt(y,B,nBTCs,zz);
    Time = toc;
    x_adct = reshape(x_hat,[H,W]);

    if CS == 2
        [x_hat,PSNR,Time] = IDA(Img,nBTCs,B,denoiser,iters,DF);        
    else
        x_hat = x_adct;
    end
    
    switch FrameType
        case 1 % Ref1
            Ref1   = x_adct;
            PFrame = x_adct;
        case 2 % Ref 2
            Ref2   = x_adct;
        case 3 % non Ref
            PFrame = x_adct;
    end

end