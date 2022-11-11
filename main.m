%% load LDH phantom data then use two segment model_A_closed for reconstruction


%% Initialization, load data 
clear all;
close all;
clc;

runID = 2;%run1 R1
          %run2 R2

          %_full retrospectively full k
load(['./run',num2str(runID),'/matlab.mat']);%

addpath('./sharedFunc');
addpath('./P2LA_Closed');
zxLoadPath;

resultFolder = ['./run',num2str(runID)];
mkdir(resultFolder);
titlePre = ['run',num2str(runID)];
 
R2_flag = 0;% 1, raw data is undersampled R=2; 
            % 0, raw data is fully sampled R=1 
ms_scale = 1000; %translate from ms to s
snr_corr_coeff = sqrt(2-pi/2);%correction for magnitude based SNR measurement
%% setup recon fdv, determine useful TR
pyr = 2; lac = 1;
tDomain = 3;time_dim = 4;
LDH_rawK = squeeze(k);%x,y,z,meta,t --> x,y,meta,t. There is only 1 slice


%% choose data segment for analyze
if runID == 1
    startTR = 21;%most voxel start to present pyr signal
    headLen = 10;%number of data point from startTR to peak lac of most voxels
    endTR = 60;
end

if runID == 2
    startTR = 21;%most voxel start to present pyr signal
    headLen = 8;%number of data point from startTR to peak lac of most voxels
    endTR = 60;
    R2_flag = 1;

end

R1_flag = not(R2_flag);

TR = 2; % constant repetition time1
FA = 20;% include as unknown, searching range [0 20]
nSlice = 1;% Keith's slice number is 1
tAxis = TR*(0:(endTR-startTR));

% fit_Mask = ones(myLen,myLen);
fdv.FA = FA;

fdv.fitvars={'kpl','T1Lac','FA','PInit','LInit'};
kpl=1;T1lac=2;fa=3;pInit=4;lInit=5;
fdv.knowns=  { 'T1Pyr','klp'};
fdv.knownvals=[ 53.5,     0    ];              
fdv.VarN = length(fdv.fitvars);
fdv.ntp = length(tAxis);% pa structure is x-y-t-z
fdv.NSeg=1;
fdv.NFlips=(fdv.ntp)*(fdv.NSeg);
fdv.TR=ones(1,fdv.NFlips)*TR;
fdv.taxis=cumsum(fdv.TR)-fdv.TR(1);
fdv.FlipAngle=fdv.FA*ones(2,fdv.NFlips);

% fdv.fitVIF = 0;%if 1, fit VIF first, then fit kpl
fdv.headLen = headLen;%how many data points are modeled with non-zero kpl
fdv.verbose=0;
VIFP_flag = 1;
    % 1--use view-shared undersampled pyr TS
    % 2--use raw data pyr TS

%% generate full k-space data
if mod((endTR - startTR),2)==0
    endTR = endTR -1;%to make sure there are even number of TRs
end

LDH_rawK = LDH_rawK(:,:,:,startTR:endTR);%truncate data from tp at startTR, before then are just noise
dataSize = size(LDH_rawK);
myLen = size(LDH_rawK,1);%in-plane 2D matrix size

pyr_k_raw =  squeeze(LDH_rawK(:,:,2,:));%x,y,t
if R2_flag
%     pyr_k_raw = pyr_k_raw.*exp(-sqrt(-1)*repmat((1:myLen/2)*1*pi/myLen,[myLen,1,size(pyr_k_raw,3)]));
end

lac_k_raw =  squeeze(LDH_rawK(:,:,1,:));%x,y,t
%% raw k-space
    PyrKmatTS = pyr_k_raw;
    LacKmatTS = lac_k_raw;
    PyrKmatTS_R1 = permute(PyrKmatTS,[2,1,3]);
    LacKmatTS_R1 = permute(LacKmatTS,[2,1,3]);
%% raw image space

%     pyr_I = fftshift(fftshift(ifft2(rot90(PyrKmatTS,-1)),1),2);
%     lac_I = fftshift(fftshift(ifft2(rot90(LacKmatTS,-1)),1),2);
    pyr_I = fftshift(fftshift(ifft2(PyrKmatTS_R1),1),2);
    lac_I = fftshift(fftshift(ifft2(LacKmatTS_R1),1),2);

%%
if R1_flag
    f2=figure(2);set(gcf,'Position',[200, 100, 400, 400]);
    imagescn(abs(pyr_I)); title('pyr R1 raw');
    export_fig(f2,[resultFolder,'/pyr_R1_raw.tif'],'-r200');
   
    f3=figure(3);set(gcf,'Position',[200, 100, 400, 400]);
    imagescn(abs(lac_I)); title('lac R1 raw');
    export_fig(f3,[resultFolder,'/lac_R1_raw.tif'],'-r200');
else
    f2=figure(2);set(gcf,'Position',[200, 100, 400, 400]);
    imagescn(abs(pyr_I)); title('pyr R2 raw');
    export_fig(f2,[resultFolder,'/pyr_R2_raw.tif'],'-r200');
   

    f3=figure(3);set(gcf,'Position',[200, 100, 400, 400]);
    imagescn(abs(lac_I)); title('lac R2 raw');
    export_fig(f2,[resultFolder,'/lac_R2_raw.tif'],'-r200');
end

if R2_flag
    pyr_k_full_undSamp = zeros(myLen,myLen,endTR-startTR+1);
    lac_k_full_undSamp = zeros(myLen,myLen,endTR-startTR+1);
    pyr_k_full_undSamp(:,1:2:end,1:2:end) =  squeeze(pyr_k_raw(:,:,1:2:end));%x,y,t
    pyr_k_full_undSamp(:,2:2:end,2:2:end) =  squeeze(pyr_k_raw(:,:,2:2:end));%x,y,t
    pyr_k_full_undSamp = permute(pyr_k_full_undSamp,[2,1,3]);
%     pyr_k_full_undSamp = pyr_k_full_undSamp.*exp(-sqrt(-1)*repmat((1:myLen)'*2*pi/myLen,[1,myLen,size(pyr_k_full_undSamp,3)]));
    %for some reason, R2 raw data has phase shift caused one voxel shift
    %along phase direction, correct it for now

    lac_k_full_undSamp(:,1:2:end,1:2:end) =  squeeze(lac_k_raw(:,:,1:2:end));%x,y,t
    lac_k_full_undSamp(:,2:2:end,2:2:end) =  squeeze(lac_k_raw(:,:,2:2:end));%x,y,t
    lac_k_full_undSamp = permute(lac_k_full_undSamp,[2,1,3]);
    
    pyr_k_full_sldWin = pyr_k_full_undSamp;
    lac_k_full_sldWin = lac_k_full_undSamp;
    pyr_k_full_sldWin(2:2:end,:,1:2:end)= pyr_k_full_sldWin(2:2:end,:,2:2:end);
    lac_k_full_sldWin(2:2:end,:,1:2:end)= lac_k_full_sldWin(2:2:end,:,2:2:end);
    pyr_k_full_sldWin(:,:,2:2:end)=[];
    lac_k_full_sldWin(:,:,2:2:end)=[];
    pyr_k_full_sldWin_auc = sum(pyr_k_full_sldWin,3);
    lac_k_full_sldWin_auc = sum(lac_k_full_sldWin,3);
    pyr_I_auc = fftshift(ifft2(pyr_k_full_sldWin_auc));
    lac_I_auc = fftshift(ifft2(lac_k_full_sldWin_auc));

    f11=figure(11); set(gcf,'Position',[200, 100, 800, 400]);

    subplot(1,2,1);imagesc(abs(pyr_I_auc)); axis image ; title('pyr sldwin auc');
    subplot(1,2,2);imagesc(abs(lac_I_auc)); axis image ; title('lac sldwin auc');
    % subplot(1,2,1);imagesc(abs(fftshift(fftshift(ifft2(pyr_k_full_sldWin_auc),2),1))); axis image ; title('pyr sldwin auc');
    % subplot(1,2,2);imagesc(abs(fftshift(fftshift(ifft2(lac_k_full_sldWin_auc),2),1))); axis image ; title('lac sldwin auc');
    sgtitle(['run',num2str(runID)]);
    export_fig(f11,[resultFolder,'/raw_slcwin_auc.tif'],'-r200'); 
    
    f2=figure(21);set(gcf,'Position',[200, 100, 400, 400]);
    pyr_I_sldwin = fftshift(fftshift(ifft2((pyr_k_full_sldWin)),1),2);
    imagescn(abs(pyr_I_sldwin));
    export_fig(f2,[resultFolder,'/pyr_slcwin_dyn.tif'],'-r200'); 
    
    f3=figure(31);set(gcf,'Position',[200, 100, 400, 400]);
    lac_I_sldwin = fftshift(fftshift(ifft2((lac_k_full_sldWin)),1),2);
    imagescn(abs(lac_I_sldwin));
%     axis image xy; 
    export_fig(f3,[resultFolder,'/lac_slcwin_dyn.tif'],'-r200'); 
else
    pyr_I_auc = sum(pyr_I,3);
    lac_I_auc = sum(lac_I,3);

    f11=figure(11); set(gcf,'Position',[200, 100, 800, 400]);
    subplot(1,2,1);imagesc(abs(pyr_I_auc)); axis image ; title('pyr auc');
    subplot(1,2,2);imagesc(abs(lac_I_auc)); axis image ; title('lac auc');
    sgtitle(['run',num2str(runID)]);
    export_fig(f11,[resultFolder,'/raw_auc.tif'],'-r200'); 
end

bkg_rot = permute(bkg,[2,1,3]);
f4=figure(4);set(gcf,'Position',[200, 100, 400, 400]);imagescn(abs(bkg_rot));axis image;
export_fig(f4,[resultFolder,'/anat.tif'],'-r200'); 

%% R2 k-space
if R2_flag
    PyrKmatTS_R = pyr_k_full_undSamp;                   
    LacKmatTS_R = lac_k_full_undSamp;  
else
    PyrKmatTS_R = PyrKmatTS_R1;                   
    LacKmatTS_R = LacKmatTS_R1;  
end
%% R2 img space
if R2_flag
    pyr_I = fftshift(fftshift(ifft2(pyr_k_full_sldWin),1),2);
    lac_I = fftshift(fftshift(ifft2(lac_k_full_sldWin),1),2);
    pyr_I_R2 = fftshift(fftshift(ifft2(pyr_k_full_undSamp),1),2);
    lac_I_R2 = fftshift(fftshift(ifft2(lac_k_full_undSamp),1),2);
   
    fit_Mask = ones([size(lac_I_sldwin,1),size(lac_I_sldwin,2)]);
else
    fit_Mask = ones([size(lac_I,1),size(lac_I,2)]);
end
%% find peak lac SNR vxl
peakSnrPyrID = find( abs(pyr_I) == max(abs(pyr_I(:))) );
[peakSnrPyrCoorX, peakSnrPyrCoorY, peakSnrPyrCoorT ]= ind2sub([size(pyr_I)],peakSnrPyrID);
peakSnrPyrTS = squeeze(pyr_I(peakSnrPyrCoorX, peakSnrPyrCoorY, :));

peakSnrLacID = find( abs(lac_I) == max(abs(lac_I(:))) );
[peakSnrLacCoorX, peakSnrLacCoorY, peakSnrLacCoorT ]= ind2sub([size(lac_I)],peakSnrLacID);
peakSnrLacTS = squeeze(lac_I(peakSnrLacCoorX, peakSnrLacCoorY, :));%Lac tSNR of the vxl with peak Pyr AUC 
peakSnrPyrAtPeakLacLocTS = squeeze(pyr_I(peakSnrLacCoorX, peakSnrLacCoorY, :));

if R2_flag
    % peakSnrPyrTS_R2 = squeeze(pyr_I_R2(peakSnrPyrCoorX, peakSnrPyrCoorY, :));
    peakSnrPyrTS_R2 = squeeze(pyr_I_R2(peakSnrLacCoorX, peakSnrLacCoorY, :));
    peakSnrLacTS_R2 = squeeze(lac_I_R2(peakSnrLacCoorX, peakSnrLacCoorY, :));
end
Lac_AUC = abs(sum(lac_I,3));
% corner_index=[1:(ceil(myLen/4)+1),1:(ceil(myLen/4)+1)];
% lac_mk = Lac_AUC>(mean(Lac_AUC(corner_index))+3*std(Lac_AUC(corner_index)));
lac_mk = Lac_AUC>(mean(Lac_AUC(:))+2*std(Lac_AUC(:)));

Pyr_AUC = abs(sum(pyr_I,3));
pyr_mk = Pyr_AUC>(mean(Pyr_AUC(:))+2*std(Pyr_AUC(:)));
recon_mk = (pyr_mk + lac_mk)>0;

%% visually exam raw data elbow point that LDH stopped working
pyr_rawK = pyr_I.*repmat(recon_mk,[1,1,size(pyr_I,3)]);
pyr_rawK_2D = reshape(pyr_rawK,[size(pyr_rawK,1)*size(pyr_rawK,2),size(pyr_rawK,3)]);
figure;plot(abs(pyr_rawK_2D)');

lac_rawK = lac_I.*repmat(recon_mk,[1,1,size(lac_I,3)]);
lac_rawK_2D = reshape(lac_rawK,[size(lac_rawK,1)*size(lac_rawK,2),size(lac_rawK,3)]);
figure;plot(abs(lac_rawK_2D)');

%%

pyr_I_mkd = pyr_I.*repmat(recon_mk,[1,1,size(pyr_I,3)]);
peakSnrPyrID = find( abs(pyr_I_mkd) == max(abs(pyr_I_mkd(:))) );
[peakSnrPyrCoorX, peakSnrPyrCoorY, peakSnrPyrCoorT ]= ind2sub([size(pyr_I_mkd)],peakSnrPyrID);
peakSnrPyrTS = squeeze(pyr_I_mkd(peakSnrPyrCoorX, peakSnrPyrCoorY, :));

lac_I_mkd = lac_I.*repmat(recon_mk,[1,1,size(lac_I,3)]);
peakSnrLacID = find( abs(lac_I_mkd) == max(abs(lac_I_mkd(:))) );
[peakSnrLacCoorX, peakSnrLacCoorY, peakSnrLacCoorT ]= ind2sub([size(lac_I_mkd)],peakSnrLacID);
peakSnrLacTS = squeeze(lac_I_mkd(peakSnrLacCoorX, peakSnrLacCoorY, :));

pyr_ns = std(abs(squeeze(pyr_I(1:3,1:3, peakSnrPyrCoorT))),[],'all');
lac_ns = std(abs(squeeze(lac_I(1:3,1:3, peakSnrLacCoorT))),[],'all');
corner_lac_vl = squeeze(lac_I(1,1, :));
peak_Pyr_vl_SNR = abs(max(peakSnrPyrTS)./std(abs(corner_lac_vl)))
peak_Lac_vl_SNR = abs(max(peakSnrLacTS)./std(abs(corner_lac_vl)))
peak_Pyr_at_Peak_Lac_vl_SNR = abs(max(peakSnrPyrAtPeakLacLocTS)./std(abs(corner_lac_vl)))

% peak_Pyr_vl_SNR = abs(max(peakSnrPyrTS))./pyr_ns*snr_corr_coeff
% peak_Lac_vl_SNR = abs(max(peakSnrLacTS))./lac_ns*snr_corr_coeff

[~,~,pyr_vl_SNR_range] = find(abs(max(pyr_I_mkd,[],3))./pyr_ns*snr_corr_coeff);
[~,~,lac_vl_SNR_range] = find(abs(max(lac_I_mkd,[],3))./lac_ns*snr_corr_coeff);

lac_mk_id = find(lac_mk);
Pyr_vl_SNR = abs(max(pyr_I_mkd,[],3))./pyr_ns*snr_corr_coeff;
Pyr_vl_SNR_array = Pyr_vl_SNR(lac_mk_id);
Pyr_vl_mean = mean(Pyr_vl_SNR_array)
Pyr_vl_std = std(Pyr_vl_SNR_array)

Pyr_vl_SNR_array_short=Pyr_vl_SNR_array(1:end-2);
mean(Pyr_vl_SNR_array_short)
std(Pyr_vl_SNR_array_short)

Lac_vl_SNR = abs(max(lac_I_mkd,[],3))./lac_ns*snr_corr_coeff;
Lac_vl_SNR_array = Lac_vl_SNR(lac_mk_id);
Lac_vl_mean = mean(Lac_vl_SNR_array)
Lac_vl_std = std(Lac_vl_SNR_array)
%%
if R2_flag
figure(31);
    subplot(2,1,1);plot(abs(peakSnrPyrTS_R2));title('peak Pyr signal');
    subplot(2,1,2);plot(abs(peakSnrLacTS_R2));title('peak Lac signal');
end

% figure;imagesc(abs(max(pyr_I,[],3)));movegui(gcf,[100,100]) % manually view peak Signal 

% fdv.VIFP = peakSnrPyrTS_R2';
fdv.VIFP = peakSnrPyrTS;
fdv.mask = fit_Mask;
fdv.mask2D = reshape(fit_Mask,[myLen.^2,1]);
fdv.data=zeros(2,fdv.NFlips,myLen*myLen);%meta,time,loc

%% Set up fdv structure:
[CMatR1, CMatR2, CMatR4, CMatR8] = undSampCMat_syn(myLen);
if R2_flag
    fdvUnd_R=fdv;
    fdvUnd_R.R=CMatR2.R;
    fdvUnd_R.CenterN=CMatR2.CenterN;
    fdvUnd_R.CenterC=CMatR2.CenterC;       
    fdvUnd_R.CMat=CMatR2.CMat;
else
    fdvUnd_R=fdv;
    fdvUnd_R.R=CMatR1.R;
    fdvUnd_R.CenterN=CMatR1.CenterN;
    fdvUnd_R.CenterC=CMatR1.CenterC;       
    fdvUnd_R.CMat=CMatR1.CMat;
end
%% initialize est result matrix
jbopts=optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','MaxFunctionEvaluations',1e10,'display','off');

LB = zeros(myLen,myLen,fdv.VarN);% fdv.fitvars={'kpl','VIFscale','ymax','a','tmax','t0'};
UB = ones(myLen,myLen,fdv.VarN);
   
%%  fdv.fitvars={'kpl','T1Lac','FA','LInit','PInit'};
LB(:,:,1) = 1e-5;
UB(:,:,1) = 1;
LB(:,:,2) = 40;
UB(:,:,2) = 53;
LB(:,:,3) = 15;
UB(:,:,3) = 20; 
LB(:,:,4) = 0;
UB(:,:,4) = 1e5;     
LB(:,:,5) = 0;
UB(:,:,5) = 5e5;    

bestresid=Inf;
fits_Mat = zeros(myLen.^2, fdv.VarN);
%% loop for every slice
for sliceID=1:nSlice

    paramInit=ones(myLen,myLen,fdv.VarN);
    LB2D = reshape(LB,[myLen.^2,fdv.VarN]);
    UB2D = reshape(UB,[myLen.^2,fdv.VarN]);
    fits_R = zeros(myLen.^2,fdv.VarN);      

    for timeID = 1:fdv.ntp
        PyrImatTS_R(:,:,timeID)=fftshift(ifft2(squeeze(PyrKmatTS_R(:,:,timeID)))); 
        LacImatTS_R(:,:,timeID)=fftshift(ifft2(squeeze(LacKmatTS_R(:,:,timeID)))); 
    end
%%  assign value to undersampled data
    fdv.mask = fit_Mask(:,:,sliceID) ;% for this testing, we just use 1 slice, but eventually will upgrade to 3D volume 
    fdvUnd_R.mask = fdv.mask;    
    fdvUnd_R.data(1,:,:) = permute(reshape(PyrImatTS_R,[myLen*myLen,fdv.ntp]),[2,1]);%meta,time,loc
    fdvUnd_R.data(2,:,:) = permute(reshape(LacImatTS_R,[myLen*myLen,fdv.ntp]),[2,1]);
    fdvUnd_R.pyrAngMap = zeros(myLen, myLen);
    fdvUnd_R.lacAngMap = zeros(myLen, myLen);        
    fdvUnd_R.pyrInitMap = zeros(myLen, myLen);     
    fdvUnd_R.lacInitMap = zeros(myLen, myLen);             
    fdvUnd_R.Len = myLen;
%     fdvUnd_R.KmatTS_pyr = pyr_k_full_undSamp;
%     fdvUnd_R.KmatTS_lac = lac_k_full_undSamp;
    fdvUnd_R.KmatTS_pyr = PyrKmatTS_R;
    fdvUnd_R.KmatTS_lac = LacKmatTS_R;
    fdvUnd_R = zxUndSmp_VIF_fit(fdvUnd_R, VIFP_flag);       
%% Init_Guess, use random
    paramInit = reshape(paramInit,[myLen.^2, fdv.VarN]);
    for parmID = 1: fdv.VarN
        tmp =  abs(randn(myLen.^2,1));
        paramInit(:,parmID)=tmp/max(tmp)*(UB2D(1,parmID)-LB2D(1,parmID))+LB2D(1,parmID);
    end

    paramInit = reshape(paramInit,[myLen,myLen,fdv.VarN]);        
    paramInit2D = reshape(paramInit,[myLen.^2,fdv.VarN]); 
%% Fit_Method

    pInit_Mat = zeros(myLen,myLen);%store the 1st point of each TS
    lInit_Mat = zeros(myLen,myLen);
    pPhase_Mat = zeros(myLen,myLen);%store the 1st point of each TS
    lPhase_Mat = zeros(myLen,myLen);
    
    %% reconstruction of 
    Guess_R = [0.05,30,15,rand,rand];
%     Guess_R1 = [0.0005,30,15,rand,rand];
    paramInit2D = repmat(Guess_R,[myLen.^2,1]);%depends on the parm_flag, cp paste the init_Guess to all voxel for R2 and R4
    vId_base = reshape(1:myLen.^2,[myLen,myLen]);
    vId_base = vId_base(1:myLen/fdvUnd_R.R,:);
    vId_array = zeros(fdvUnd_R.R,1);
    for tmpId = 1:myLen.^2/fdvUnd_R.R
        vId_array(1)=vId_base(tmpId);
        for IDtmp = 1:fdvUnd_R.R-1
            vId_array(1+IDtmp) = vId_array(1) + myLen*IDtmp/fdvUnd_R.R;
        end

        fdvFit_R = fdvUnd_R;

        mask_array = fdv.mask(vId_array);
        mask_array_nonZInd = find(mask_array);%only fit the non-zero voxels
        mask_array_ZInd = find(mask_array==0);%indentify bg zero voxels

        GuessGen = paramInit2D(vId_array(mask_array_nonZInd),:) ;
        LBGen =  [LB2D(vId_array(mask_array_nonZInd),:) repmat(-pi, length(mask_array_nonZInd), 2)];       
        UBGen =  [UB2D(vId_array(mask_array_nonZInd),:) repmat(pi, length(mask_array_nonZInd), 2)];                     

        fdvFit_R.mask_array_nonZInd = mask_array_nonZInd;%such as 1,2 out of 1:4
        fdvFit_R.mask_array_ZInd = mask_array_ZInd;% such as 3,4 out of 1:4
        fdvFit_R.voxFitInd = vId_array(mask_array_nonZInd);         

        if sum(mask_array_nonZInd) > 0
            pyr_tmp = squeeze(fdvUnd_R.data(1,:,vId_array(fdvUnd_R.CenterC)));
            lac_tmp = squeeze(fdvUnd_R.data(2,:,vId_array(fdvUnd_R.CenterC)));

            fdvFit_R.data = ([pyr_tmp;lac_tmp]); %determine location of the voxel in the undersampled image   
            %% use fitted phase and init value 
            vxlPyrInit = abs(zxUnaliasTsInit_R_general(pyr_tmp,fdvUnd_R.CMat));
            vxlPyrInitMked = vxlPyrInit.*recon_mk(vId_array);
            fdvFit_R.pyrInit = vxlPyrInitMked;

            vxlLacInit = abs(zxUnaliasTsInit_R_general(lac_tmp,fdvUnd_R.CMat));
            vxlLacInitMked = vxlLacInit.*recon_mk(vId_array);
            fdvFit_R.lacInit = vxlLacInitMked;
            fdvFit_R.recon_mk = recon_mk(vId_array);

            fdvUnd_R.pyrInitMap(vId_array(mask_array_nonZInd)) = fdvFit_R.pyrInit(mask_array_nonZInd);
            fdvUnd_R.lacInitMap(vId_array(mask_array_nonZInd)) = fdvFit_R.lacInit(mask_array_nonZInd);

            %% use the ground truth from non-undersampled data
            % 'kpl','T1Lac','FA','LInit','PInit'
            GuessGen(:,end-1) = abs(fdvFit_R.pyrInit(mask_array_nonZInd));%assign PInit                   
            GuessGen(:,end)   = abs(fdvFit_R.lacInit(mask_array_nonZInd));%assign LInit 
            GuessGen = [GuessGen zeros(length(mask_array_nonZInd), 2)];%assign phase for pyr & lac
            [fitGenNonZInd,resid] = lsqnonlin(@(x) P2LAClosed_2Seg_Err_R_general_IntegPhase(x,fdvFit_R),GuessGen, LBGen, UBGen, jbopts);       
            fitGenNonZInd(:,end-1:end) = [];%delete pyr and lac phase           
                         
            fits_R(vId_array(mask_array_nonZInd),:)=fitGenNonZInd;
            fits_R(vId_array(mask_array_ZInd),1)=0;                   
        else
            fits_R(vId_array,1) = 0;
        end     
    end     


fits_Mat(:,:,sliceID) = fits_R;  
end

%% display pyr/lac TS for publication
disp_TR = length(tAxis); %display portion TS, start from T0, ignore later noise 
% mask_Mat_pyr = repmat(fit_Mask,[1,1,disp_TR]);
% mask_Mat_lac = repmat(lac_mk,[1,1,disp_TR]);

f_width =600;
numCol=4;
interval=5;
fontSize=12;

[R_Pyr_fit_TS, R_Lac_fit_TS] = zxCompose_PL_TS_2Seg(fits_R,fdvUnd_R,disp_TR);
R_PL_fit_TS = zeros(myLen,myLen*2,disp_TR);
R_PL_fit_TS(:,1:myLen,:)=R_Pyr_fit_TS;
R_PL_fit_TS(:,myLen+1:myLen*2,:)=R_Lac_fit_TS*3;

f321=figure(321);set(gcf,'color','w','Position',[200, 100, 600,330]);
PyrMax = max(abs(R_Pyr_fit_TS(:)));
LacMax = max(abs(R_Lac_fit_TS(:)));
for col = 1:numCol
    subtightplot(2,numCol,col);         imagesc(abs(R_Pyr_fit_TS(:,:,((col-1)*interval+1))),[0,PyrMax]);axis off image;colormap gray;title(['t=',num2str(((col-1)*interval+1)*fdv.TR(1)),'s'],'FontSize', fontSize);
    subtightplot(2,numCol,col+numCol);  imagesc(abs(R_Lac_fit_TS(:,:,((col-1)*interval+1))),[0,LacMax]);axis off image;colormap gray;        
end
sgtitle('Estimated');   
export_fig(f321,[resultFolder,'/Recon_TS_4figs.tif'],'-r200');

    %%

    f320=figure(320);set(gcf,'color','w','Position',[200, 100, f_width,330]);
    % imagescn(abs(R1PL_fit_TS)*scaleRatio);
    PyrMax = max(abs(PyrImatTS_R(:)));
    LacMax = max(abs(LacImatTS_R(:)));
    for col = 1:numCol
        subtightplot(2,numCol,col);         imagesc(abs(PyrImatTS_R(:,:,((col-1)*interval+1))),[0,PyrMax]);axis off image;colormap gray;title(['t=',num2str(((col-1)*interval+1)*fdv.TR(1)),'s'],'FontSize', fontSize);
        subtightplot(2,numCol,col+numCol);  imagesc(abs(LacImatTS_R(:,:,((col-1)*interval+1))),[0,LacMax]);axis off image;colormap gray;
        
    end
    sgtitle('Undersampled');
    export_fig(f320,[resultFolder,'/UndSamp_TS_4figs.tif'],'-r200');    
    
    f3201=figure(3201);set(gcf,'color','w','Position',[200, 100, f_width,330]);
    pyr_R2_fits_2D = reshape(R_Pyr_fit_TS,[myLen*myLen,length(tAxis)]);
    lac_R2_fits_2D = reshape(R_Lac_fit_TS,[myLen*myLen,length(tAxis)]);
    subplot(2,1,1);plot(pyr_R2_fits_2D');title('R2 pyr fit dyn');
    subplot(2,1,2);plot(lac_R2_fits_2D');title('R2 lac fit dyn');
    export_fig(f3201,[resultFolder,'/fit_dyn.tif'],'-r200');   
    
    f3202=figure(3202);set(gcf,'color','w','Position',[200, 100, f_width,330]);
    % imagescn(abs(R1PL_fit_TS)*scaleRatio);
    PyrMax = max(abs(pyr_I(:)));
    LacMax = max(abs(lac_I(:)));
    for col = 1:numCol
        subtightplot(2,numCol,col);         imagesc(abs(pyr_I(:,:,((col-1)*interval+1))),[0,PyrMax]);axis off image;colormap gray;title(['t=',num2str(((col-1)*interval+1)*fdv.TR(1)),'s'],'FontSize', fontSize);
        subtightplot(2,numCol,col+numCol);  imagesc(abs(lac_I(:,:,((col-1)*interval+1))),[0,LacMax]);axis off image;colormap gray;       
    end
    sgtitle('R=2 sliwin dyn');
    export_fig(f3202,[resultFolder,'/sliwin_dyn.tif'],'-r200');  
    
%% reshape 2D mat(x*y,t) to 3D mat(x,y,t)
fitResult_R = reshape(fits_Mat,[myLen,myLen,fdv.VarN,nSlice]);

f22=figure(22);set(gcf,'Position',[200, 100, 400, 400]);
subtightplot(2,3,1);imagesc(squeeze(fitResult_R(:,:,kpl)));colorbar;title('kpl');axis image;
subtightplot(2,3,3);imagesc(squeeze(fitResult_R(:,:,T1lac)));colorbar;title('T1Lac');axis image;
subtightplot(2,3,4);imagesc(squeeze(fitResult_R(:,:,fa)));colorbar;title('FA');axis image;
subtightplot(2,3,5);imagesc(squeeze(abs(fitResult_R(:,:,pInit))));colorbar;title('P0');axis image;
subtightplot(2,3,6);imagesc(squeeze(abs(fitResult_R(:,:,lInit))));colorbar;title('L0');axis image;
subtightplot(2,3,2);imagesc(squeeze(fitResult_R(:,:,kpl)).*recon_mk);colorbar;title('kpl within recon mask');axis image;axis off;

% f22=figure(22);set(gcf,'Position',[200, 100, 600, 200]);
% subtightplot(1,3,1);imagesc(squeeze(fitResult_R2(:,:,kpl)));colorbar;title('kpl');axis image;
% subtightplot(1,3,2);imagesc(squeeze(abs(fitResult_R2(:,:,pInit))));colorbar;title('P0');axis image;
% subtightplot(1,3,3);imagesc(squeeze(abs(fitResult_R2(:,:,lInit))));colorbar;title('L0');axis image;

fitKpl_R = squeeze(fitResult_R(:,:,kpl)).*recon_mk;
% fitKplR2 = fitKplR2(:);
% fitKpl_R = fitKpl_R(find(fitKpl_R(:)>(1e-4)))
sgtitle('fit results');
export_fig(f22,[resultFolder,'/fit_result.tif'],'-transparent','-r200');   
dlmwrite([resultFolder,'/kpl_mk.txt'], fitKpl_R);
%% save result
save([resultFolder,'/',titlePre]);
