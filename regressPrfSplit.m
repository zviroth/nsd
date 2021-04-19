%uses file from prfSampleModel_nsd2.m

%ignore first and last pyramid levels.

close all
clear all
tic
bandpass = 1; bandMin = 3; bandMax = 6;

isub=1;

toZscore=0;

nsplits=2;
prffolder = '~/NSD/prfsample/';
betasfolder = ['~/NSD/sub' num2str(isub) '_betas_func1pt8mm/'];
% stimfilename = fullfile(folder,'nsdsynthetic_colorstimuli_subj01.hdf5');
nsdfolder = '~/NSD/';
roifolder = ['~/NSD/sub' num2str(isub) '_betas_func1pt8mm/'];
visualRoisFile = fullfile(roifolder,'prf-visualrois.nii');%V1v, V1d, V2v, V2d, V3v, V3d, and hV4
visRoiData = niftiread(visualRoisFile);
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4'};

visRoiData = visRoiData(:);

% for toZscore=0:0
%     for nsdOrSynth=0:1
%         if nsdOrSynth==1
load(fullfile(prffolder,['prfSampleStim_v1_sub' num2str(isub) '.mat']),'prfSampleLevOri','prfSampleLev',...
    'rois','allImgs','numLevels','numOrientations','interpImgSize','backgroundSize','pixPerDeg',...
    'roiPrf');

if bandpass
%     for iroi=rois
%         prfSampleLevOri{iroi} = prfSampleLevOri{iroi}(:,:,2:numLevels-1,:);
%         prfSampleLev{iroi} = prfSampleLev{iroi}(:,:,2:numLevels-1);
%     end
%     numLevels = numLevels-2;
    for iroi=rois
        prfSampleLevOri{iroi} = prfSampleLevOri{iroi}(:,:,bandMin:bandMax,:);
        prfSampleLev{iroi} = prfSampleLev{iroi}(:,:,bandMin:bandMax);
    end
    numLevels = bandMax-bandMin+1;
end

for iroi=rois; roiBetas{iroi}=[]; end
for isession=1:40
    betasfilename = fullfile(betasfolder,['betas_session' num2str(isession,'%02.f') '.nii']);
    betas = niftiread(betasfilename);
    betas = cast(betas,'double');
    betas = betas/300;
    betas=reshape(betas,[],size(betas,4));
    for iroi=rois
        if toZscore
            roiBetas{iroi} = [roiBetas{iroi} zscore(betas(visRoiData==iroi,:),0,2)];
        else
            roiBetas{iroi} = [roiBetas{iroi} betas(visRoiData==iroi,:)];
        end
    end
end


nsdDesignFilename = fullfile(nsdfolder, 'nsd_expdesign.mat');
nsdDesign = load(nsdDesignFilename);
nsdDesign.masterordering;%for each of 30000 trials, what is corresponding image (out of 10000 images)
subDesign = nsdDesign.subjectim(isub,nsdDesign.masterordering);%for each of 30000 trials, what is corresponding image (out of 73000 images)
[imgTrials, imgNum] = ismember(subDesign, allImgs);%logical array

imgTrialSum = cumsum(imgTrials);
totalTrials = max(imgTrialSum);%3000
midTrial = find(imgTrialSum==ceil(totalTrials/2));
splitImgTrials = repmat(imgTrials,2,1);
splitImgTrials(1,midTrial+1:end) = zeros;
splitImgTrials(2,1:midTrial) = zeros;

r2 = cell(length(rois),1);
r2ori = cell(length(rois),1);
r2oriSplit = cell(length(rois),1);
r2split = cell(length(rois),1);
%         else
%             load(fullfile(prffolder,['prfSampleSynth_v1_sub' num2str(isub) '.mat']),'prfSampleLevOri','prfSampleLev',...
%                 'rois','allImgs','numLevels','numOrientations','interpImgSize','backgroundSize','pixPerDeg',...
%                 'roiPrf');
%             synthfolder = ['~/NSD/sub' num2str(isub) '_synth_func1pt8mm/'];
%             betasfilename = fullfile(synthfolder,'betas_nsdsynthetic.nii');
%             betas = niftiread(betasfilename);
%             betas = cast(betas,'double');
%             betas = betas/300;
%
%             synthDesignFilename = fullfile(nsdfolder, 'nsdsynthetic_expdesign.mat');
%             synthDesign = load(synthDesignFilename);
%             synthDesign.masterordering;
%
%             [imgTrials, imgNum] = ismember(synthDesign.masterordering, allImgs);%logical array
%             betas=reshape(betas,[],size(betas,4));
%
%             for iroi=rois
%                 roiBetas{iroi} = betas(visRoiData==iroi,:);
%             end
%         end


for iroi=rois
    %     roiBetas{iroi} = betas(visRoiData==iroi,:);
    
    nvox(iroi) = size(roiBetas{iroi},1);
    voxOriResidual{iroi} = zeros(nsplits, nvox(iroi),ceil(totalTrials/2));
    voxResidual{iroi} = zeros(nsplits, nvox(iroi),ceil(totalTrials/2));
    voxOriResidualSplit{iroi} = zeros(nsplits, nvox(iroi),ceil(totalTrials/2));
    voxResidualSplit{iroi} = zeros(nsplits, nvox(iroi),ceil(totalTrials/2));
    voxOriCoef{iroi} = zeros(nsplits, nvox(iroi),numLevels*numOrientations+1);
    voxCoef{iroi} = zeros(nsplits, nvox(iroi),numLevels+1);
    voxPredOriCoef{iroi} = zeros(nsplits, nvox(iroi),numLevels*numOrientations+1);
    voxPredCoef{iroi} = zeros(nsplits, nvox(iroi),numLevels+1);
    voxOriPredOriCoef{iroi} = zeros(nsplits, nvox(iroi),numLevels*numOrientations+1);
    %get model coefficients for each voxel, within each split
    for isplit=1:nsplits
        imgTrials = splitImgTrials(isplit,:);
        for ivox=1:nvox(iroi)
            voxBetas = roiBetas{iroi}(ivox,imgTrials>0)';
            voxPrfSample = squeeze(prfSampleLev{iroi}(imgNum(imgTrials>0),ivox,:));
            %add constant predictor
            voxPrfSample(:,end+1) = ones;
            %voxePrfSample*coeff = voxBetas;
            voxCoef{iroi}(isplit,ivox,:) = voxPrfSample\voxBetas;%check vox 144 in first ROI
            
            voxPrfOriSample = squeeze(prfSampleLevOri{iroi}(imgNum(imgTrials>0),ivox,:,:));
            voxPrfOriSample = reshape(voxPrfOriSample,[],numLevels*numOrientations);
            %add constant predictor
            voxPrfOriSample(:,end+1) = ones;
            voxOriCoef{iroi}(isplit,ivox,:) = voxPrfOriSample\voxBetas;%check vox 144 in first ROI
        end
    end
    
    %compute R^2 for both models, within splits, and between splits
    for isplit=1:nsplits
        imgTrials = splitImgTrials(isplit,:);
        for ivox=1:nvox(iroi)
            voxBetas = roiBetas{iroi}(ivox,imgTrials>0)';
            voxPrfSample = squeeze(prfSampleLev{iroi}(imgNum(imgTrials>0),ivox,:));
            voxPrfSample(:,end+1) = ones;
            voxPrfOriSample = squeeze(prfSampleLevOri{iroi}(imgNum(imgTrials>0),ivox,:,:));
            voxPrfOriSample = reshape(voxPrfOriSample,[],numLevels*numOrientations);
            voxPrfOriSample(:,end+1) = ones;
            
            voxOriResidual{iroi}(isplit,ivox,:) = voxBetas' - squeeze(voxOriCoef{iroi}(isplit,ivox,:))'*voxPrfOriSample';
            voxResidual{iroi}(isplit,ivox,:) = voxBetas' - squeeze(voxCoef{iroi}(isplit,ivox,:))'*voxPrfSample';
            
            voxOriResidualSplit{iroi}(isplit,ivox,:) = voxBetas' - squeeze(voxOriCoef{iroi}(nsplits-isplit+1,ivox,:))'*voxPrfOriSample';
            voxResidualSplit{iroi}(isplit,ivox,:) = voxBetas' - squeeze(voxCoef{iroi}(nsplits-isplit+1,ivox,:))'*voxPrfSample';
            
            %regress vignetting predicted timecourse with orientation model
            voxPred = squeeze(voxCoef{iroi}(isplit,ivox,:))'*voxPrfSample';
            voxOriPred = squeeze(voxOriCoef{iroi}(isplit,ivox,:))'*voxPrfOriSample';
            %             voxPrfOriSample = squeeze(prfSampleLevOri{iroi}(imgNum(imgTrials>0),ivox,:,:));
            %             voxPrfOriSample = reshape(voxPrfOriSample,[],numLevels*numOrientations);
            %             voxPrfOriSample(:,end+1) = ones;%add constant predictor
            voxPredOriCoef{iroi}(isplit,ivox,:) = voxPrfOriSample\voxPred';
            voxOriPredOriCoef{iroi}(isplit,ivox,:) = voxPrfOriSample\voxOriPred';
            
            
            %             [B,BINT,R,RINT,STATS] = regress(voxBetas,voxPrfSample);
            
        end
        %r2 within split
        r2{iroi}(isplit,:) = 1 - (rssq(voxResidual{iroi}(isplit,:,:),3)'./rssq(roiBetas{iroi}(:,imgTrials>0),2)).^2;
        r2ori{iroi}(isplit,:) = 1 - (rssq(voxOriResidual{iroi}(isplit,:,:),3)'./rssq(roiBetas{iroi}(:,imgTrials>0),2)).^2;
        
        %r2 between splits
        r2split{iroi}(isplit,:) = 1 - (rssq(voxResidualSplit{iroi}(isplit,:,:),3)'./rssq(roiBetas{iroi}(:,imgTrials>0),2)).^2;
        r2oriSplit{iroi}(isplit,:) = 1 - (rssq(voxOriResidualSplit{iroi}(isplit,:,:),3)'./rssq(roiBetas{iroi}(:,imgTrials>0),2)).^2;
        
    end
end
% voxResidual{iroi}(isplit,ivox,:) = voxBetas' - squeeze(voxCoef{iroi}(isplit,ivox,:))'*voxPrfSample';

%%
figure(1); clf
rows=2; cols=5;
isubplot=0;

isubplot = isubplot+1; subplot(rows,cols,isubplot)
temp = reshape(prfSampleLev{iroi},size(prfSampleLev{iroi},1)*size(prfSampleLev{iroi},2),size(prfSampleLev{iroi},3));
plot(mean(temp));
title('mean power in pRF'); xlabel('pyramid level');

%
isubplot = isubplot+1; subplot(rows,cols,isubplot)
%

for isplit=1:nsplits
    isubplot = isubplot+1; subplot(rows,cols,isubplot)
    plot(squeeze(voxCoef{iroi}(isplit,:,:))')
    title('mean voxel coefficients'); xlabel('pyramid level');
end
%
for isplit=1:nsplits
    isubplot = isubplot+1; subplot(rows,cols,isubplot)
    plot(squeeze(voxCoef{iroi}(isplit,:,:)))
    title('mean voxel coefficients'); xlabel('voxel #');
end


isubplot = isubplot+1; subplot(rows,cols,isubplot)
hist(r2{iroi}(:)); hold on; vline(median(r2{iroi}(:)));
title('r^2, vignetting model');



isubplot = isubplot+1; subplot(rows,cols,isubplot)
hist(r2ori{iroi}(:)); hold on; vline(median(r2ori{iroi}(:)));
title('r^2, full model');



isubplot = isubplot+1; subplot(rows,cols,isubplot)
hist(r2split{iroi}(:)); hold on; vline(median(r2split{iroi}(:)));
title('r^2 SPLIT, vignetting model');

isubplot = isubplot+1; subplot(rows,cols,isubplot)
hist(r2oriSplit{iroi}(:)); hold on; vline(median(r2oriSplit{iroi}(:)));
title('r^2 SPLIT, full model');

% if nsdOrSynth==1b
nsd.voxResidual = voxResidual;
nsd.voxOriResidual = voxOriResidual;
nsd.voxResidualSplit = voxResidualSplit;
nsd.voxOriResidualSplit = voxOriResidualSplit;
nsd.r2 = r2;
nsd.r2ori = r2ori;
nsd.r2split = r2split;
nsd.r2oriSplit = r2oriSplit;

nsd.imgTrials = imgTrials;
nsd.imgNum = imgNum;
nsd.splitImgTrials = splitImgTrials;
nsd.midTrial = midTrial;

nsd.voxCoef = voxCoef;
nsd.voxOriCoef = voxOriCoef;

nsd.voxPredOriCoef = voxPredOriCoef;
nsd.voxOriPredOriCoef = voxOriPredOriCoef;
%         else
%             synth.roiBetas = roiBetas;
%             synth.voxOriCoef = voxOriCoef;
%             synth.voxCoef = voxCoef;
%             synth.prfSampleLevOri = prfSampleLevOri;
%             synth.prfSampleLev = prfSampleLev;
%             synth.r2 = r2;
%             synth.r2ori = r2ori;
%             synth.voxResidual = voxResidual;
%             synth.voxOriResidual = voxOriResidual;
%             synth.imgTrials = imgTrials;
%             synth.imgNum = imgNum;
% end
%     end
zscoreStr='';
if toZscore
    zscoreStr = '_zscore';
end
bandpassStr = '';
if bandpass
    bandpassStr = ['_bandpass' num2str(bandMin) 'to' num2str(bandMax)];
end
save(fullfile(prffolder,['regressPrfSplit' bandpassStr '_v1_sub' num2str(isub) zscoreStr '.mat']), 'nsd', ...
    'numLevels', 'numOrientations','rois','nvox','roiPrf','nsplits');

toc
