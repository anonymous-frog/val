% -------------------------------------------------------------------------
% VAL-VFI-VidSet6
% -------------------------------------------------------------------------

warning('OFF'); clc; clearvars; rng(10); % close all;

clear Codec_Type;

addpath(genpath("DnCNN"));
addpath matconvnet-1.0-beta25/matlab

n_DnCNN_layers = 17; %Options 17 & 20 (17 is better)
LoadNetworkWeights(n_DnCNN_layers);

fprintf('%s\n\n',datestr(now));

%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------

ADAPT_ref    = 3;       % 1,2,3 = ZZ,BBV,DD
ADAPT_nonref = 3;       % 1,2,3 = ZZ,BBV,DD
CS_ref       = 0;       % 0,1,2 = No CS, DAMP-D, IDA
CS_nonref    = 0;       % 0,1,2 = No CS, DAMP-D, IDA
Key_offset   = 0;       % MUST BE 0
Dynamic_Ref  = 0;       % 0 or 1 - Has no role to play with Video Frame Interpolation
DPCM         = 1;       % 0 or 1
VFI          = 1;       % 0 or 1
CORRECTION   = 1;       % 0 or 1
Corr_Thresh  = 25;      % 0 => prefer predicted frame  >0 => prefer non predicted non key frame

B_key        = 16;      % Block size Key
B            = 16;      % Block size non Key
DF_K         = 2;       % Damping factor for IDA Key frame
DF_nK        = 2;       % Damping factor for IDA non Key frame
Iters_K      = 20;      % Number of iterations for IDA Key frame
Iters_nK     = 20;      % Number of iterations for IDA non Key frame

Ref_Type     = 2;       % (1) THB-DD (2) THI-DD
Cur_Type     = 2;       % (1) THB-DD (2) THI-DD (3) RTDD (4) PTDD (5) MTDD

GOP_size     = 4;       % Number of frames in GOP
num_frames   = 16;      % Number of frames in SEQ

denoiser     = 'DnCNN'; % Denoiser 'DnCNN'

DnCNNF        = 0;       % 0 => no filter >1 => filter before VFI, filter on Rx with DnCNN value

DISPLAY      = 2;       % 0,1,2 no display, display all, display 2 frames
WAIT         = 0;       % If DISPLAY == 2, WAIT == 1 => click on image
TIMING       = 0;       % 0,1 no timing info, show timing info
MONTAGE      = 0;       % show all frames in one image
SAVE         = 0;       % save results on file

key_subrate  = 0.5;     % Key subrate
% non Key subrate so that overall GOP subrate is same as 0.7 (key) and 0.1 (non key)
%s ubrate = 0.1+(0.7-key_subrate)/(GOP_size-1);  
subrate      = ((GOP_size*0.175)-key_subrate)/(GOP_size-1);

ADAPT        = [0,1,9]; % ZZ, BBV, DD
ADAPTS       = {'ZZ','BBV','DD'};

disp([ADAPT_ref ADAPT_nonref CS_ref CS_nonref Key_offset Dynamic_Ref DPCM VFI]);

n_DnCNN_layers = 17; %Options 17 & 20 (17 is better)
LoadNetworkWeights(n_DnCNN_layers);

DIR = 'Data/Videos/';

Videos = {'paris_cif', 'foreman_cif', 'coastguard_cif', ...
    'hall_monitor_cif', 'mobile_cif', 'news_cif'};     

Set = [1 2 3 4 5 6];

PSNRA = []; SSIMA = []; TIMEA = [];
TIME1 = 0;  TIME3 = 0;

RES = zeros(numel(Set),...  % SETP = Sequence Number
    num_frames,...          % n Number of frames
    3 ...                   % PSNR,SSIM,TIME PSNR,SSIM,TIME
    );

GOP_Seq = 1:GOP_size;

[~,zz] = Zz(B);
t_mcp=0; t_mcp_tcs=0; t_dpcm=0; t_corr=0;

SETP = 0;
for SET = Set
    
    SEQ = cell(1,num_frames);
    
    SETP = SETP + 1;
    PSNR = zeros(1,num_frames);
    SSIM = zeros(1,num_frames);
    TIME = 0;
    
    sequence_name = [DIR Videos{SET} '.y4m'];
    f = fopen(sequence_name);
    
    %fprintf('-----------------\n');
    fprintf('%s \n',Videos{SET});
    %fprintf('-----------------\n');
    
    for k = 1:GOP_size:num_frames
               
        if SET == 5
            frames = readframes(f, GOP_size+1, k, 10); % Avoid error in mobile_cif (SET == 5)
        else
            frames = readframes(f, GOP_size+1, k, 0);
        end
        
        %frames = readframes_sif(f, GOP_size, k); % SIF for stefan_sif.y4m
        [num_rows, num_cols] = size(frames{1});
        
        n_store = 0; % use to skip already read Key frames in MID Key mode
        for n = GOP_Seq
            
            Time_Start = clock;
            
            nfp = k + n - 1;
            
            if (n == GOP_Seq(1)) || (n == (GOP_Seq(1)+GOP_size))
                
                if n ~= n_store
                    n_store = n;
                    %---------------------------------- ARRef1 -----------------
                    reference_frame = frames{n};
                    [H, W] = size(reference_frame);
                    Ref1 = reference_frame;
                    
                    CFi = 1/key_subrate;
                    
                    [ARRef1, nBTCs, T1, CFa, PSNRs] = Codec_Type( ...
                        reference_frame, 1, Ref_Type, ...
                        B_key, ADAPT(ADAPT_ref), CFi, CS_ref, ...
                        Iters_K, DF_K, denoiser);
                    
                    TIME = TIME + T1;
                                        
                    % No time expended since these would have been received
                    TCs_Key1 = zeros(H,W);
                    for r = 1:B:H
                        for c = 1:B:W
                            Patch = ARRef1(r:r+B-1,c:c+B-1);
                            TCs_Key1(r:r+B-1,c:c+B-1)=dct2(Patch);
                        end
                    end
                    
                    Psnr = psnr(ARRef1,reference_frame,255);
                    PSNR(nfp) = Psnr;
                    Ssim = msssim(ARRef1,reference_frame);
                    SSIM(nfp) = Ssim;
                    
                    fprintf('%6.4f %5.2f %6.4f %6.4f',1/CFa, Psnr, Ssim, T1);
                    
                    T3 = 0;
                    if TIMING == 1
                        fprintf('T1=%6.4f T3=%6.4f \n',T1,T3);
                    else
                        fprintf('\n');
                    end
                    
                    if DISPLAY == 1
                        subplot(231); imshow(ARRef1,[0 255]);
                        title(['Ref ' num2str(n+k-1) ' ' num2str(Psnr)]);
                        drawnow;
                    end
                    
                    if DISPLAY == 2
                        
                        figure(1);
                        
                        subplot(131);
                        imshow(imresize(reference_frame,2),[0 255]);
                        Title = sprintf('Key frame %02d', ...
                            n+k-1);
                        title(Title);
                        
                        subplot(132);
                        imshow(nBTCs(1:B:H,1:B:W),[]);                        
                        
                        subplot(133);
                        imshow(imresize(ARRef1,2),[0 255]);
                        RES(SETP,nfp,1) = Psnr;
                        RES(SETP,nfp,2) = Ssim;
                        Title = sprintf('Coded Key frame %02d  %5.2f dB %6.4f', ...
                            n+k-1,Psnr, Ssim);
                        title(Title);
                        
                        Seq = Videos{SET}; Seqp = strfind(Seq,'_'); Seq(Seqp) = '-';
                        if WAIT == 0
                            SGTitle = sprintf('%s sequence',Seq);
                        else
                            SGTitle = sprintf('%s sequence.\n Click on Image to proceed.', ...
                                Seq);
                        end
                        drawnow;
                        if WAIT == 1
                            waitforbuttonpress;
                        end
                        
                        SEQ{1,nfp} = ARRef1;
                        
                    end
                    
                    %---------------------------------- ARRef2 -----------------
                    reference_frame = frames{n+GOP_size};
                    [H, W] = size(reference_frame);

                    [ARRef2, ~ , T1, CFa, PSNRs] = Codec_Type( ... 
                        reference_frame, 2, Ref_Type, ...
                        B_key, ADAPT(ADAPT_ref), CFi, CS_ref, ...
                        Iters_K, DF_K, denoiser);
                    
                    % No time expended since these would have been received
                    TCs_Key1 = zeros(H,W);
                    for r = 1:B:H
                        for c = 1:B:W
                            Patch = ARRef2(r:r+B-1,c:c+B-1);
                            TCs_Key1(r:r+B-1,c:c+B-1)=dct2(Patch);
                        end
                    end
                    
                end
                                
                tic
                
                ARRef1c(:,:,1) = ARRef1;
                ARRef1c(:,:,2) = ARRef1;
                ARRef1c(:,:,3) = ARRef1;
                
                ARRef2c(:,:,1) = ARRef2;
                ARRef2c(:,:,2) = ARRef2;
                ARRef2c(:,:,3) = ARRef2;
                
                imwrite(uint8(ARRef1c),'DAIN-master/VAL/Key/Key_1.png');
                imwrite(uint8(ARRef2c),'DAIN-master/VAL/Key/Key_2.png');
                
                cd 'DAIN-master'
                [status, result] = ...
                    system(['python -W ignore DAIN_vfi.py --netName DAIN_slowmotion --time_step ' ...
                    num2str(1/GOP_size)]);
                cd '..'
                
                t_mcp = toc;
                               
            end
            
            if (n ~= 1 + Key_offset) && (n ~= 9 + Key_offset)
                
                %---------------------------------- ARCur -----------------
                current_frame = frames{n};
                [H, W] = size(current_frame);
                Cur = current_frame;
                
                CFi = 1/subrate;
                
                if Cur_Type == 4    
                    if (n == 2)
                        C_Type = 2; % THI
                    else
                        C_Type = 4; % Pframe-TDD
                    end
                else
                    C_Type = Cur_Type;
                end
                                                              
                [ARCur, nBTCs, T1, CFa, PSNRs] = Codec_Type( ...
                    current_frame, 3, C_Type, ...                   
                    B, ADAPT(ADAPT_nonref), CFi, CS_nonref, ...
                    Iters_nK, DF_nK, denoiser);
                                               
                TIME1 = TIME1 + T1;
                                
                
                if DISPLAY == 1
                    subplot(234); imshow(ARCur,[0 255]);
                    title('ABCS Current'); drawnow;
                end
                
                % Calculate the TCs of the BBV or DD recon Cur Frame
                % No time expended since these would have been received
                TCs_Cur = zeros(H,W);
                for r = 1:B:H
                    for c = 1:B:W
                        Patch = ARCur(r:r+B-1,c:c+B-1);
                        D = dct2(Patch);
                        D(zz(nBTCs(r,c)+1:end))=0;
                        TCs_Cur(r:r+B-1,c:c+B-1)=D;
                    end
                end
                
                %==========================================================
                %                      VFI
                %==========================================================                
                                               
                if VFI == 1
                    fname = sprintf('%04d.png',n-1);
                    VFI_Ref = imread(['DAIN-master/VAL/nonKey/' fname]);
                    full_name = ['Data/Non_Key/' Videos{SET} '/' ...
                        Videos{SET} '_' num2str(nfp) '.png'];
                    imwrite(VFI_Ref,full_name);
                    VFI_Ref = double(rgb2gray(VFI_Ref));
                else
                    VFI_Ref = (ARRef1+ARRef2)/2;
                end
                
                if DISPLAY == 1
                    subplot(232), imshow(VFI_Ref,[0 255]);
                    title('VFI Ref'); drawnow;
                end
                
                % Calculate the TCs of the Motion Compensated Ref Frame
                tic
                TCs_MC_Ref = zeros(H,W);
                for r = 1:B:H
                    for c = 1:B:W
                        Patch = VFI_Ref(r:r+B-1,c:c+B-1);
                        TCs_MC_Ref(r:r+B-1,c:c+B-1)=dct2(Patch);
                    end
                end
                t_mcp_tcs = toc;
                
                %==========================================================
                %                       DPCM
                %==========================================================
                
                if DPCM == 1
                    tic
                    Mask = (abs(TCs_Cur)==0);
                    TCs = TCs_MC_Ref.*Mask + TCs_Cur;
                    
                    Rx = zeros(H,W);
                    for r = 1:B:H
                        for c = 1:B:W
                            D = TCs(r:r+B-1,c:c+B-1);
                            Rx(r:r+B-1,c:c+B-1)=idct2(D);
                        end
                    end
                    t_dpcm = toc;
                else
                    Rx = VFI_Ref;
                end
                
                if (DPCM == 0) && (VFI == 0)
                    Rx = ARCur;
                end
                
                Psnr = psnr(Rx,current_frame,255);
                Ssim = msssim(Rx,current_frame);
                
                if DISPLAY == 1
                    subplot(235); imshow(Rx,[0 255]);
                    title(['VFI-DPCM ' num2str(n+k-1) ' ' num2str(Psnr)]);
                end
                
                
                % =========================================================
                %                   Correction
                % =========================================================
                
                tic
                
                if CORRECTION == 1
                    
                    Err = Rx-VFI_Ref;
                    Mask = abs(Err)>Corr_Thresh;
                    Rx_Corr = Rx.*(1-Mask) + ARCur.*Mask;
                    
                else
                    
                    Rx_Corr = Rx;
                    
                end
                
                %----------------------------------------------------------
                %                   FILTER
                %----------------------------------------------------------
                              
                if DnCNNF > 1
                    
                    Y = Rx_Corr(:);
                    if SET == 8
                        X = denoise(Y,sqrt(100),288,352,'DnCNN');
                    else
                        X = denoise(Y,sqrt(BM3DF),288,352,'DnCNN');
                    end
                    Rx_Corr = reshape(X,[288 352]);
                    
                end
                
                t_corr = toc;
                
                %----------------------------------------------------------
                %                   STATS and DISPLAY
                %----------------------------------------------------------
                
                if DISPLAY == 1
                    subplot(233); imshow(Mask);
                end
                
                Psnr = psnr(Rx_Corr,current_frame,255);  PSNR(nfp) = Psnr;
                Ssim = msssim(Rx_Corr,current_frame);    SSIM(nfp) = Ssim;
                
                if DISPLAY == 1
                    subplot(236); imshow(Rx_Corr,[0 255]);
                    title(['Corrected ' num2str(n+k-1) ' ' num2str(Psnr)]);
                    drawnow; % waitforbuttonpress;
                end
                
                if DISPLAY == 2
                    
                    figure(1);
                    
                    subplot(131);
                    imshow(imresize(current_frame,2),[0 255]);
                    Title = sprintf('Cur frame %02d',n+k-1);
                    title(Title);
                    
                    subplot(132);
                    imshow(nBTCs(1:B:H,1:B:W),[]);
                    
                    subplot(133);
                    imshow(imresize(Rx_Corr,2),[0 255]);
                    RES(SETP,nfp,1) = Psnr;
                    RES(SETP,nfp,2) = Ssim;
                    Title = sprintf('Coded Cur frame %02d  %5.2f dB %6.4f', ...
                        n+k-1,Psnr, Ssim);
                    title(Title);
                    
                    Seq = Videos{SET}; Seqp = strfind(Seq,'_'); Seq(Seqp) = '-';
                    if WAIT == 0
                        SGTitle = sprintf('%s sequence',Seq);
                    else
                        SGTitle = sprintf('%s sequence.\n Click on Image to proceed.', ...
                            Seq);
                    end
                    
                    drawnow;
                    if WAIT == 1
                        waitforbuttonpress;
                    end
                    
                    SEQ{1,nfp} = Rx_Corr;
                    
                end
                
                
                if Dynamic_Ref == 1
                    ARRef1 = (VFI_Ref+ARRef1)/2;
                    ARRef2 = (VFI_Ref+ARRef2)/2;
                end
                
                time = t_mcp/(GOP_size - 1) + t_mcp_tcs + t_dpcm + t_corr;
                
                fprintf('%6.4f %5.2f %6.4f ',1/CFa, Psnr, Ssim);
                
                if TIMING == 1
                    fprintf('T1=%6.4f T3=%6.4f t_mcp=%6.4f ',...
                        T1,T3,t_mcp);
                    fprintf('t_mcp_tcs=%6.4f t_dpcm=%6.4f t_corr=%6.4f\n',...
                        t_mcp_tcs,t_dpcm,t_corr);
                else
                    fprintf('%6.4f\n',T1 + time);
                end
                
                
                TIME1 = TIME1 + time;
                RES(SETP,nfp,3) = T1 + time;
                
            end
            
            Time = etime(clock, Time_Start);
            TIME = TIME + Time;
            
        end
        
    end
    
    fprintf('------------\n');
    fprintf('%5.2f %6.4f %8.2fs\n\n',mean(PSNR),mean(SSIM),TIME);
    
    PSNRA = [PSNRA mean(PSNR)];
    SSIMA = [SSIMA mean(SSIM)];
    TIMEA = [TIMEA TIME];
    
    
    if MONTAGE == 1
        figure; montage(SEQ,'DisplayRange',[0 255]);
    end
    
    if SAVE == 1
        saveas(gcf,[Videos{SET} '.png']);
    end
    
end

disp([ADAPT_ref ADAPT_nonref CS_ref CS_nonref Key_offset Dynamic_Ref DPCM VFI]);

disp([PSNRA' SSIMA']);
fprintf('   -----------------\n');
fprintf('   %5.2f      %6.4f\n',mean(PSNRA),mean(SSIMA));

mTime = mean(TIMEA);

Encoding_Time_per_frame    = mTime/16;
Decoding_Time_per_frame_D1 = TIME1/num_frames/6;
Decoding_Time_per_frame_RT = TIME3/num_frames/6;

fprintf('\nEncoding Time per Frame = %7.4f\n',Encoding_Time_per_frame);
fprintf('Decoding Time per Frame = %7.4f for Denoiser 1\n',Decoding_Time_per_frame_D1);
fprintf('Decoding Time per Frame = %7.4f for Real Time \n\n',Decoding_Time_per_frame_RT);

fprintf('%s\n\n',datestr(now));



