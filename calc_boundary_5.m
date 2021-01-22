function [nBTCs,CS,nCS] = jz_calc_boundary_5(Img,B,CF)
%
% INPUTS
%   Img       Image
%   B         Block size
%   CF        Compression Factor N/M > 1
%
% OUTPUTS
%   nBTCs     number of TCs per block stored in nBTCs(1:B:H-1,1:B:W-1)
%   CS        actual contour pixels sensed
%   nCS       number of pixel differences sensed sum(sum(CS~=0))/2

% If CF>32 then we collect no BBV since Step and Start become Inf

    ns = min(B/2,fix(B/CF));
    
    [H,W] = size(Img);
    NB = H*W/B^2; % Number of Blocks
    NS = ns*2*NB + ns*fix((H-1)/B) + ns*fix((W-1)/B); % Expected BBV measurements
    
    Img(H+1,1:W) = Img(H,1:W);
    Img(1:H,W+1) = Img(1:H,W);
    Img(H+1,W+1) = Img(H,W);
    
    H = H + 1;
    W = W + 1;
    
    Step  = max(round(B/ns),2);
  
    Start = floor(Step/2);
    nb    = 2*ns;
    
    M = zeros(H,W);     % Measurements
    F = zeros(H,W);     % Filter with pixel Positions
    BBV = zeros(H,W);   % Block Boundary Variation
    nBTCs = zeros(H,W); % Number of TCs in each Block
    
    % Capture initial TV Measurments
    for R = 1:B:H
        for C = 1:B:W
            
            % Collect one block row
            if C == W
                % Do not collect a block row which would extend off image
            else
                for c = C+Start-1:Step:C+B-1
                    try
                        M(R,c) = abs(Img(R,c)-Img(R,c+1));
                    catch
                        fprintf('Error!\n');
                    end
                    F(R,c) = 1;
                    F(R,c+1) = 1;
                end
            end
            
            % Collect one block column
            if R == H
                % Do not collect a block column  which would extend off image
            else
                for r = R+Start-1:Step:R+B-1
                    M(r,C) = abs(Img(r,C)-Img(r+1,C));
                    F(r,C) = 1;
                    F(r+1,C) = 1;
                end
            end
            
        end
    end
    
    % Calculate Block Boundary Variation
    for R = 1:B:H-1
        for C = 1:B:W-1
            BBV(R,C) = BBV(R,C) + sum(M(R:R+B,C));
            BBV(R,C) = BBV(R,C) + sum(M(R:R+B,C+B));
            BBV(R,C) = BBV(R,C) + sum(M(R,C:C+B));
            BBV(R,C) = BBV(R,C) + sum(M(R+B,C:C+B));
            BBV(R,C) = BBV(R,C) + 1; % Avoids BBV = 0
        end
    end
    
    % Calculate number of TCs in each block = nBTCs[]
    
    TBBV    = sum(BBV(:));  % Total Block Boundary Variation
    TTCs    = fix((H-1)*(W-1)/CF) - sum(F(:))/2; % Total Number of TCs in image MINUS BBV measurements
    %TTCs   = fix((H-1)*(W-1)/CF) - NS; % Total Number of TCs in image MINUS BBV measurements
    ETCs    = 0; % Extra TCs in total in blocks with TCs > B^2
    nBETCs  = 0; % Number of blocks with TCs > B^2
    cTCs    = 0; % Current number of TCs
    
    for R = 1:B:H-1
        for C = 1:B:W-1
            nBTCs(R,C) = fix(TTCs * BBV(R,C)/TBBV);
            if nBTCs(R,C)>B^2
                nBETCs = nBETCs + 1;
                ETCs = ETCs + nBTCs(R,C) - B^2;
                nBTCs(R,C) = B^2;     
            end
            cTCs = cTCs + nBTCs(R,C);
        end
    end
        
    diff_per_block = round((TTCs - cTCs)/NB);
       
    for R = 1:B:H-1
        for C = 1:B:W-1
            if nBTCs(R,C)<B^2
                nBTCs(R,C) = nBTCs(R,C) + diff_per_block;
                if nBTCs(R,C)>B^2
                    nBTCs(R,C) = B^2;
                end
            end
        end
    end
    
    CS = Img.*F; CS(:,W-1)=CS(:,W); CS(H-1,:) = CS(H,:);
    CS = CS(1:H-1,1:W-1);
    nCS = sum(CS(:)~=0);
    

end