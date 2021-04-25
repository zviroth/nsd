%uses files from nsdStim.m

close all
clear all

isub=1;

visualRegion = 4;%V1,V2,V3,V4

% allImgs = 1:1000;%1:1000;
nsdfolder = '~/NSD/';
nsdDesignFilename = fullfile(nsdfolder, 'nsd_expdesign.mat');
nsdDesign = load(nsdDesignFilename);
% allImgs = nsdDesign.sharedix;
allImgs = nsdDesign.subjectim(isub,nsdDesign.masterordering);%indices of all 10000 images used for this subject
allImgs = unique(allImgs);


colors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250], [0.4940, 0.1840, 0.5560], ...
    [0.4660, 0.6740, 0.1880],[0.3010, 0.7450, 0.9330], [0.6350, 0.0780, 0.1840]	};

interpImgSize = 714;
backgroundSize = 1024;
imgScaling = 0.5;

bandwidth = 1;
numOrientations = 8;
dims = [backgroundSize backgroundSize];
dims = dims*imgScaling;
numLevels = maxLevel(dims,bandwidth);

switch visualRegion
    case 1
        rois=1:2;%1:2;
    case 2
        rois=3:4;
    case 3
        rois=5:6;
    case 4
        rois=7;
end

pixPerDeg = imgScaling*714/8.4;%=85
degPerPix = 8.4/(714*imgScaling);
x = -backgroundSize/2+0.5:backgroundSize/2-0.5;%correct?
y = -backgroundSize/2+0.5:backgroundSize/2-0.5;
[X,Y] = meshgrid(x,y);
% pyramidfolder = '~/NSD/stimuli/pyramid/';
pyramidfolder = ['/Volumes/nsd_sub' num2str(isub) '/pyramid/'];
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4'};



% nrois=7;

betasfolder = ['~/NSD/sub' num2str(isub) '_betas_func1pt8mm/'];

angFile = fullfile(betasfolder,'prf_angle.nii');
eccFile = fullfile(betasfolder,'prf_eccentricity.nii');
sizeFile = fullfile(betasfolder,'prf_size.nii');
gainFile = fullfile(betasfolder,'prf_gain.nii');
expFile = fullfile(betasfolder,'prf_exponent.nii');
r2file = fullfile(betasfolder,'prf_R2.nii');

% angInfo = niftiinfo(angFile);
angData = niftiread(angFile);
eccData = niftiread(eccFile);
sizeData = niftiread(sizeFile);
gainData = niftiread(gainFile);
expData = niftiread(expFile);
r2Data = niftiread(r2file);

visualRoisFile = fullfile(betasfolder,'prf-visualrois.nii');%V1v, V1d, V2v, V2d, V3v, V3d, and hV4
visRoiData = niftiread(visualRoisFile);


% profile on
tic



%%

roiData = visRoiData;

roiPrf = cell(length(rois),1);
prfSampleLev = cell(length(rois),1);
prfSampleLevOri = cell(length(rois),1);


for roinum=1:length(rois)
    iroi = rois(roinum);
    ecc = eccData(roiData(:)==iroi);
    ang = angData(roiData(:)==iroi);
    sz = sizeData(roiData(:)==iroi);% The definition of pRF size is one standard deviation of a Gaussian that describes the response of the model to point stimuli. Note that this definition of pRF size takes into account the nonlinear summation behavior of the pRF model and is not the same as the ?sigma? parameter used in the pRF model.
    exponent = expData(roiData(:)==iroi);
    gain = gainData(roiData(:)==iroi);
    nvox = length(ecc);
    r2 = r2Data(roiData(:)==iroi);
    
    roiPrf{roinum}.ecc = ecc;
    roiPrf{roinum}.ang = ang;
    roiPrf{roinum}.sz = sz;
    roiPrf{roinum}.exponent = exponent;
    roiPrf{roinum}.gain = gain;
    roiPrf{roinum}.r2 = r2;
    
    sz = sz*pixPerDeg;
    [x0, y0] = pol2cart(ang*2*pi/360,ecc*pixPerDeg);
    
    roiPrf{roinum}.x = x0;
    roiPrf{roinum}.y = y0;
    
%     prfSampleLev{iroi} = zeros(nvox,length(allImgs),numLevels);
%     prfSampleLevOri{iroi} = zeros(nvox,length(allImgs),numLevels, numOrientations);
    
    prfSampleLevRoi = zeros(length(allImgs),nvox,numLevels);
    prfSampleLevOriRoi = zeros(length(allImgs),nvox,numLevels, numOrientations);
    
    %loop through synthetic images
    parfor iimg=1:length(allImgs)
        ['roi: ' num2str(iroi) ', image: ' num2str(iimg)]
        imgNum = allImgs(iimg);
        pyramidfilename = ['pyrImg' num2str(imgNum) '.mat'];
        data = load(fullfile(pyramidfolder, pyramidfilename),'sumOri','modelOri');
%         load(fullfile(pyramidfolder, pyramidfilename), 'interpImgSize','backgroundSize','numOrientations',...
%             'bandwidth','dims','bigImg','sumOri','modelOri','numLevels');
        imgPrfSampleLev = zeros(nvox,numLevels);
        imgPrfSampleLevOri = zeros(nvox,numLevels,numOrientations);
        %loop through voxels
        for ivox=1:nvox
%         for ivox=1:nvox
            %sample for each pyramid level (summed across orientations!)
            voxPrfSampleLev = zeros(numLevels,1);
            voxPrfSampleLevOri = zeros(numLevels,numOrientations);
            
            G = exp(-((X-x0(ivox)).^2+(Y-y0(ivox)).^2)/(2*sz(ivox)^2));
            
            for ilev = 1:numLevels
                voxPrfSampleOri = zeros(numOrientations,1);
%                 voxPrfSampleLev(ilev) = prfResp(sumOri{ilev},X,Y,x0(ivox),y0(ivox),sz(ivox),exponent(ivox),gain(ivox));
                voxPrfSampleLev(ilev) = gain(ivox)*((data.sumOri{ilev}(:)'*G(:))^exponent(ivox));
                %                 prfSampleLev{iroi}(ivox,iimg,ilev) = prfResp(sumOri{ilev},X,Y,x0(ivox),y0(ivox),sz(ivox),exponent(ivox),gain(ivox));
                for iori=1:numOrientations
%                     voxPrfSampleOri(iori) = prfResp(modelOri{ilev}(iori,:,:),X,Y,x0(ivox),y0(ivox),sz(ivox),exponent(ivox),gain(ivox));
                    temp = data.modelOri{ilev}(iori,:,:);
                    voxPrfSampleOri(iori) = gain(ivox)*((temp(:)'*G(:))^exponent(ivox));
                    %                     prfSampleLevOri{iroi}(ivox,iimg,ilev,iori) = prfResp(modelOri{ilev}(iori,:,:),X,Y,x0(ivox),y0(ivox),sz(ivox),exponent(ivox),gain(ivox));
                end
                voxPrfSampleLevOri(ilev,:) = voxPrfSampleOri;
            end
            imgPrfSampleLevOri(ivox,:,:) = voxPrfSampleLevOri;
            imgPrfSampleLev(ivox,:) = voxPrfSampleLev;
        end
        prfSampleLevRoi(iimg,:,:) = imgPrfSampleLev;
        prfSampleLevOriRoi(iimg,:,:,:) = imgPrfSampleLevOri;
    end
    prfSampleLev{roinum} = prfSampleLevRoi;
    prfSampleLevOri{roinum} = prfSampleLevOriRoi;
end

% end
prffolder = '~/NSD/prfsample/';
save(fullfile(prffolder,['prfSampleStim_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']),'prfSampleLevOri','prfSampleLev',...
    'rois','allImgs','numLevels','numOrientations','interpImgSize','backgroundSize','pixPerDeg',...
    'roiPrf','-v7.3');

%save it also on the external drive
backupfolder = ['/Volumes/nsd_sub' num2str(isub) '/prfSample/'];
save(fullfile(backupfolder,['prfSampleStim_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']),'prfSampleLevOri','prfSampleLev',...
    'rois','allImgs','numLevels','numOrientations','interpImgSize','backgroundSize','pixPerDeg',...
    'roiPrf','-v7.3');

toc

% profile viewer
%%
% plot example voxels
figure
exampleVox = [535 456 290];
rows=length(exampleVox)+1;
cols=1+numLevels;
%originals
subplot(rows,cols,1)
pyramidfilename = ['pyrImg' num2str(allImgs(end)) '.mat'];
load(fullfile(pyramidfolder, pyramidfilename),'bigImg','sumOri');
imagesc(bigImg);
for ilev=1:numLevels
    subplot(rows,cols,1+ilev)
    imagesc(sumOri{ilev});
end
%example voxels
for i=1:length(exampleVox)
    ivox=exampleVox(i);
    G = exp(-((X-x0(ivox)).^2+(Y-y0(ivox)).^2)/(2*sz(ivox)^2));
    %original image
    SxG = bigImg.*G;
    subplot(rows,cols,cols+(i-1)*cols+1)
    imagesc(SxG);
    %pyramid levels
    for ilev=1:numLevels
        SxG = sumOri{ilev}.*G;
        subplot(rows,cols,cols+(i-1)*cols+1+ilev)
        imagesc(SxG);
    end
end
for isubplot=1:rows*cols
    subplot(rows,cols,isubplot)
    axis image
    colormap gray
    axis off
end
set(gcf,'position',[100 140 1350 620]);

%% CSS pRF from Kay et al., J Neurophysiology, 2013
function r = prfResp(S,x,y,x0,y0,sigma, n, g)
G = exp(-((x-x0).^2+(y-y0).^2)/(2*sigma^2));
% r = g*(sum(S(:).*G(:))^n);
r = g*((S(:)'*G(:))^n);
end