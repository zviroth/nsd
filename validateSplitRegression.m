%uses data from regressPrfSplit.m

close all
clear all

isub=1;


bandpass = 1; bandMin = 1; bandMax = 6;

toZscore=0;
nperms=100;
prffolder = '~/NSD/prfsample/';
zscoreStr='';
if toZscore
    zscoreStr = '_zscore';
end

absAvg = 0;

bandpassStr = '';
if bandpass
    bandpassStr = ['_bandpass' num2str(bandMin) 'to' num2str(bandMax)];
end
load(fullfile(prffolder,['regressPrfSplit' bandpassStr '_v1_sub' num2str(isub) zscoreStr '.mat']), 'nsd', ...
    'numLevels', 'numOrientations','rois','nvox','roiPrf','nsplits');
% nsplits=2;
interpImgSize = 714;
backgroundSize = 1024;
%%
figure(1); clf
rows=length(rois); cols=6;
isubplot=0;

for iroi=rois
    
        isubplot = isubplot+1; subplot(rows,cols,isubplot)
    hist(nsd.r2ori{iroi}(:)); hold on; vline(median(nsd.r2ori{iroi}(:)));
    title('within, full model');
    isubplot = isubplot+1; subplot(rows,cols,isubplot)
    hist(nsd.r2{iroi}(:)); hold on; vline(median(nsd.r2{iroi}(:)));
    title('within, vignetting model');
    isubplot = isubplot+1; subplot(rows,cols,isubplot)
    temp = nsd.r2ori{iroi}(:) - nsd.r2{iroi}(:);
    hist(temp(:)); hold on; vline(median(temp(:)));
    title('within, full-vignetting');
    
    
        isubplot = isubplot+1; subplot(rows,cols,isubplot)
    hist(nsd.r2oriSplit{iroi}(:)); hold on; vline(median(nsd.r2oriSplit{iroi}(:)));
    title('split, full model');
    isubplot = isubplot+1; subplot(rows,cols,isubplot)
    hist(nsd.r2split{iroi}(:)); hold on; vline(median(nsd.r2split{iroi}(:)));
    title('split, vignetting model');
    
    isubplot = isubplot+1; subplot(rows,cols,isubplot)
    temp = nsd.r2oriSplit{iroi}(:) - nsd.r2split{iroi}(:);
    hist(temp(:)); hold on; vline(median(temp(:)));
    title('split, full-vignetting');
end
for isubplot=1:rows*cols
   subplot(rows,cols,isubplot)
   xlabel('R^2');
   ylabel('#voxels');
end
    
figure(2); clf;
rows=length(rois); cols=4;
isubplot=0;

for iroi=rois
    coefCorr = corr(squeeze(nsd.voxCoef{iroi}(1,:,1:end-1))', squeeze(nsd.voxCoef{iroi}(2,:,1:end-1))');
    
    isubplot = isubplot+1; subplot(rows,cols,isubplot)
    hist(diag(coefCorr)); hold on; vline(median(diag(coefCorr)));
    title('vignetting coef corr');
    
    isubplot = isubplot+1; subplot(rows,cols,isubplot)
    o = eye(size(coefCorr));
    permCoefCorr = coefCorr(o==0);
    hist(permCoefCorr); hold on; vline(median(permCoefCorr));
    title('permuted coef corr');
    
    coefCorrOri = corr(squeeze(nsd.voxOriCoef{iroi}(1,:,1:end-1))', squeeze(nsd.voxOriCoef{iroi}(2,:,1:end-1))');
    
    isubplot = isubplot+1; subplot(rows,cols,isubplot)
    hist(diag(coefCorrOri)); hold on; vline(median(diag(coefCorrOri)));
    title('full coef corr');
    
    isubplot = isubplot+1; subplot(rows,cols,isubplot)
    o = eye(size(coefCorrOri));
    permCoefCorrOri = coefCorrOri(o==0);%permuted across voxels
    hist(permCoefCorrOri); hold on; vline(median(permCoefCorrOri));
    title('permuted coef corr');
end
for isubplot=1:rows*cols
   subplot(rows,cols,isubplot)
   xlabel('correlation (r)');
   ylabel('#voxels');
end


%% PLOT PREFERRED ORIENTATION
figure(3); clf;
rows=length(rois); cols=4;
isubplot=0;

for iroi=rois
    for isplit=1:nsplits
        isubplot = isubplot+1; subplot(rows,cols,isubplot)
        plotOriLines(squeeze(nsd.voxOriCoef{iroi}(isplit,:,1:end-1)), roiPrf{iroi},numLevels,numOrientations)
%         plotOriLines(squeeze(nsd.voxCoef{iroi}(1,:,1:end-1)), roiPrf{iroi},numOrientations,linelength,linewidth)
        title(['split ' num2str(isplit)]);
    end
    %average across levels
    for isplit=1:nsplits
        isubplot = isubplot+1; subplot(rows,cols,isubplot)
        coef = reshape(squeeze(nsd.voxOriCoef{iroi}(isplit,:,1:end-1)),nvox(iroi),numLevels,numOrientations);
        if absAvg
            oriMeanLev = squeeze(mean(abs(coef),2));
        else
            oriMeanLev = squeeze(mean(coef,2));
        end
        plotOriLines(oriMeanLev, roiPrf{iroi},numLevels,numOrientations)
        %         plotOriLines(squeeze(nsd.voxCoef{iroi}(1,:,1:end-1)), roiPrf{iroi},numOrientations,linelength,linewidth)
        title(['split ' num2str(isplit) ', avg levels']);
    end
    
    
end

set(gcf,'position',[50 80 1300 700]);

%% PLOT PREFERRED LEVEL
figure(4); clf;
rows=length(rois); cols=6;
isubplot=0;
% maxwidth=0.02;
% maxlength = 1.2;
maxSize = 0.05;
for iroi=rois
    for isplit=1:nsplits
        isubplot = isubplot+1; subplot(rows,cols,isubplot)
        plotLevelDots(squeeze(nsd.voxOriCoef{iroi}(isplit,:,1:end-1)), roiPrf{iroi},numLevels,numOrientations)
        title(['split ' num2str(isplit) ', full']);
    end
    
    
    %vignetting model
    for isplit=1:nsplits
        isubplot = isubplot+1; subplot(rows,cols,isubplot)
        plotLevelDots(squeeze(nsd.voxCoef{iroi}(isplit,:,1:end-1)), roiPrf{iroi},numLevels,numOrientations)
        title(['split ' num2str(isplit) ', vignetting']);
    end
    
        %average across orientations
    for isplit=1:nsplits
        isubplot = isubplot+1; subplot(rows,cols,isubplot)
        coef = reshape(squeeze(nsd.voxOriCoef{iroi}(isplit,:,1:end-1)),nvox(iroi),numLevels,numOrientations);
        levMeanOri = squeeze(mean(coef,3));
        if absAvg
            levMeanOri = squeeze(mean(abs(coef),2));
        else
            levMeanOri = squeeze(mean(coef,3));
        end

        plotLevelDots(levMeanOri, roiPrf{iroi},numLevels,numOrientations)
        title(['split ' num2str(isplit) ', avg ori']);
    end
    
end

set(gcf,'position',[50 80 1300 700]);


%scale oriented scatter plot with the improvement in R2 for oriented model
%over vignetting model
temp = nsd.r2oriSplit{iroi}(:) - nsd.r2split{iroi}(:);

%% PLOT PREFERRED ORIENTATION - from vignetting data
figure(5); clf;
rows=length(rois); cols=6;
isubplot=0;

for iroi=rois
    
    %Use vignetting model coefficients to create predicted betas (by
    %definition contains no orientation information per se) and regress
    %using full orientation model
    for isplit=1:nsplits
        isubplot = isubplot+1; subplot(rows,cols,isubplot)
        plotOriLines(squeeze(nsd.voxPredOriCoef{iroi}(isplit,:,1:end-1)), roiPrf{iroi},numLevels,numOrientations)
        title(['fromVigToOri, split ' num2str(isplit)]);
    end
    %Use full orientation model coefficients to create predicted betas and regress
    %using full orientation model. Should return the original coefficients.
    for isplit=1:nsplits
        isubplot = isubplot+1; subplot(rows,cols,isubplot)
        plotOriLines(squeeze(nsd.voxOriPredOriCoef{iroi}(isplit,:,1:end-1)), roiPrf{iroi},numLevels,numOrientations)
        title(['fromOritoOri, split ' num2str(isplit)]);
    end
    
    %average across levels
    for isplit=1:nsplits
        isubplot = isubplot+1; subplot(rows,cols,isubplot)
        coef = reshape(squeeze(nsd.voxPredOriCoef{iroi}(isplit,:,1:end-1)),nvox(iroi),numLevels,numOrientations);
        oriMeanLev = squeeze(mean(coef,2));
        if absAvg
            oriMeanLev = squeeze(mean(abs(coef),2));
        else
            oriMeanLev = squeeze(mean(coef,2));
        end

        plotOriLines(oriMeanLev, roiPrf{iroi},numLevels,numOrientations)
        %         plotOriLines(squeeze(nsd.voxCoef{iroi}(1,:,1:end-1)), roiPrf{iroi},numOrientations,linelength,linewidth)
        title(['fromVigToOri, split ' num2str(isplit) ', avg levels']);
    end
end

set(gcf,'position',[50 80 1300 700]);


%% SCATTER PLOT MAX FILTER LEVEL VS ECCENTRICITY
figure(6); clf;
rows=length(rois); cols=2*nsplits;
isubplot=0;
maxwidth=0.02;
maxlength = 1.2;
maxSize = 0.4;
for iroi=rois
    r2 = roiPrf{iroi}.r2;
    dotSize = maxSize*max(0.0001,r2);
    for isplit=1:nsplits
        [tempOriMaxValue prefOriFilter] = max(squeeze(nsd.voxOriCoef{iroi}(isplit,:,1:end-1)),[],2);
        prefLevel = mod(prefOriFilter-1,numLevels)+1;
         isubplot = isubplot+1; subplot(rows,cols,isubplot);
         scatter(roiPrf{iroi}.ecc, prefLevel, dotSize, prefLevel); hold on;
         for ilev=1:numLevels
            scatter(median( roiPrf{iroi}.ecc(prefLevel==ilev)), ilev, maxSize*2*100, 'k');
         end
         caxis([1 numLevels]);
         title(['full, split ' num2str(isplit)]);
         xlabel('pRF eccentricity');
    end
    
    for isplit=1:nsplits    
        [tempOriMaxValue prefFilter] = max(squeeze(nsd.voxCoef{iroi}(isplit,:,1:end-1)),[],2);
%         prefAngle = (prefFilter-1)*(pi)/numOrientations;
         isubplot = isubplot+1; subplot(rows,cols,isubplot);
        scatter(roiPrf{iroi}.ecc, prefFilter, dotSize, prefFilter); hold on
        for ilev=1:numLevels
            scatter(median( roiPrf{iroi}.ecc(prefFilter==ilev)), ilev, maxSize*3*100, 'k');
        end
        caxis([1 numLevels]);
        title(['vignetting, split ' num2str(isplit)]);
        xlabel('pRF eccentricity');
    end
end 
set(gcf,'position',[100 100 900 600]);


%% mean preferred model level as function of eccentricity
figure(8); clf;
rows=length(rois); cols=nsplits;
isubplot = 0;
nbins = 10;
binBorders = logspace(-0.5,1,nbins+1);
binCoef = zeros(nsplits,nbins,numLevels);
for iroi=rois
    for isplit=1:nsplits
       for ibin=1:nbins
           binVoxels = roiPrf{iroi}.ecc>binBorders(ibin) & roiPrf{iroi}.ecc<=binBorders(ibin+1);
           binCoef(isplit,ibin,:) = squeeze(mean(nsd.voxCoef{iroi}(isplit,binVoxels,1:end-1),2));
       end
       isubplot = isubplot+1; subplot(rows,cols,isubplot);
       plot(squeeze(max(binCoef(isplit,:,:),[],3)));
       xlabel('eccentricity bin');
       ylabel('max pyramid level');
       title(['vignetting, split ' num2str(isplit)]);
    end
end

%% prf size as function of eccentricity
figure(9); clf;
rows=length(rois); cols=nsplits;
isubplot = 0;
nbins = 10;
binBorders = logspace(-0.5,1,nbins+1);
% binCoef = zeros(nsplits,nbins,numLevels);
for iroi=rois
    for isplit=1:nsplits
       for ibin=1:nbins
           binVoxels = roiPrf{iroi}.ecc>binBorders(ibin) & roiPrf{iroi}.ecc<=binBorders(ibin+1);
           binSz(isplit,ibin,:) = squeeze(mean(roiPrf{iroi}.sz(binVoxels)));
%            binCoef(isplit,ibin,:) = squeeze(mean(nsd.voxCoef{iroi}(isplit,binVoxels,1:end-1),2));
       end
       isubplot = isubplot+1; subplot(rows,cols,isubplot);
       plot(squeeze(mean(binSz(isplit,:,:),3)));
       xlabel('eccentricity bin');
       ylabel('mean pRf size');
    end
end

%% SCATTER PLOT MAX FILTER LEVEL VS PRF SIZE
figure(10); clf;
rows=length(rois); cols=2*nsplits;
isubplot=0;
maxwidth=0.02;
maxlength = 1.2;
maxSize = 0.4;
for iroi=rois
    r2 = roiPrf{iroi}.r2;
    dotSize = maxSize*max(0.0001,r2);
    for isplit=1:nsplits
        [tempOriMaxValue prefOriFilter] = max(squeeze(nsd.voxOriCoef{iroi}(isplit,:,1:end-1)),[],2);
        prefLevel = mod(prefOriFilter-1,numLevels)+1;
         isubplot = isubplot+1; subplot(rows,cols,isubplot);
         scatter(roiPrf{iroi}.sz, prefLevel, dotSize, prefLevel); hold on;
         for ilev=1:numLevels
            scatter(median( roiPrf{iroi}.sz(prefLevel==ilev)), ilev, maxSize*2*100, 'k');
         end
         caxis([1 numLevels]);
         title(['full, split ' num2str(isplit)]);
         xlabel('pRF size');
    end
    
    for isplit=1:nsplits    
        [tempOriMaxValue prefFilter] = max(squeeze(nsd.voxCoef{iroi}(isplit,:,1:end-1)),[],2);
%         prefAngle = (prefFilter-1)*(pi)/numOrientations;
         isubplot = isubplot+1; subplot(rows,cols,isubplot);
        scatter(roiPrf{iroi}.sz, prefFilter, dotSize, prefFilter); hold on
        for ilev=1:numLevels
            scatter(median( roiPrf{iroi}.sz(prefFilter==ilev)), ilev, maxSize*3*100, 'k');
        end
        caxis([1 numLevels]);
        title(['vignetting, split ' num2str(isplit)]);
        xlabel('pRF size');
    end
end 
set(gcf,'position',[100 100 900 600]);


%%
figure;
iroi=1;
isplit=1;
plotLevelDots(squeeze(nsd.voxCoef{iroi}(isplit,:,1:end-1)), roiPrf{iroi},numLevels,numOrientations)
title(['roi ' num2str(iroi) ' split ' num2str(isplit)]);
set(gcf,'position',[300 300 500 500]);



%%
%for testing the reshape.
[ori lev] = meshgrid(1:numOrientations, 1:numLevels); 
ordervec = 1:32;
reshape(ori,[],numLevels*numOrientations);
ceil(ordervec/numLevels); %preferred orientation
reshape(lev,[],numLevels*numOrientations);
mod((ordervec-1),numLevels)+1; %preferred level



%%
function res = plotOriLines(oriCoef, roiPrf, numLevels,numOrientations)

numvox = size(oriCoef,1);

maxwidth=0.02;
maxlength = 1.2;
lineWidth = maxwidth*max(0.0001,roiPrf.r2);
lineLength = maxlength*max(0.0001,roiPrf.r2);



cMap = parula(256);
% [tempOriMaxValue prefFilter] = max(abs(oriCoef),[],2);
[tempOriMaxValue prefFilter] = max(oriCoef,[],2);
if size(oriCoef,2)>numOrientations
%     prefOri = mod((prefFilter-1),numOrientations)+1;
    prefOri = ceil(prefFilter/numLevels);
    prefAngle = (prefOri-1)*(pi)/numOrientations;
else
%     [tempOriMaxValue prefOri] = max(abs(oriCoef),[],2);
     prefAngle = (prefFilter-1)*(pi)/numOrientations;
end

for ivox=1:numvox
    drawOriLine(roiPrf.x(ivox), roiPrf.y(ivox), prefAngle(ivox), lineLength(ivox), lineWidth(ivox), cMap(1+floor(prefAngle(ivox)*256/pi),:))
    hold on
end

interpImgSize = 714;
backgroundSize = 1024;
set(gca,'XTick',[]);
set(gca,'YTick',[]);
linecolor = [0.5 0.5 0.5];
line([0 0], [-backgroundSize backgroundSize],'color',linecolor);
line([-backgroundSize backgroundSize],[0 0], 'color',linecolor);
line([-interpImgSize -interpImgSize], [-interpImgSize interpImgSize], 'color',linecolor);
line([interpImgSize interpImgSize], [-interpImgSize interpImgSize], 'color',linecolor);
line([-interpImgSize interpImgSize],[interpImgSize interpImgSize], 'color',linecolor);
line([-interpImgSize interpImgSize],[-interpImgSize -interpImgSize], 'color',linecolor);
xlim([-interpImgSize interpImgSize]); ylim([-interpImgSize interpImgSize]);


axis square

end

%%
function res = plotLevelDots(oriCoef, roiPrf, numLevels, numOrientations)

numvox = size(oriCoef,1);


cMap = parula(256);
[tempOriMaxValue, prefFilter] = max(oriCoef,[],2);

if size(oriCoef,2)>numLevels
%     prefLevel = ceil(prefFilter/numOrientations);
    prefLevel = mod(prefFilter-1,numLevels)+1;
else
    prefLevel = prefFilter;
end




% for ivox=1:numvox
%     h=plot(roiPrf.x(ivox), roiPrf.y(ivox),'o','MarkerSize',dotSize(ivox),'MarkerEdgeColor',cMap(floor(prefLevel(ivox)*256/numLevels),:));
% %     h=plot(roiPrf.x(ivox), roiPrf.y(ivox),'o','MarkerSize',dotSize(ivox),'MarkerEdgeColor',cMap(floor(prefLevel(ivox)*256/numLevels),:), 'MarkerFaceColor',cMap(floor(prefLevel(ivox)*256/numLevels),:));
% %     h=plot(roiPrf.x(ivox), roiPrf.y(ivox),'o','MarkerSize',100,'MarkerEdgeColor','k', 'MarkerFaceColor','k');
% %     drawOriLine(roiPrf.x(ivox), roiPrf.y(ivox), 'o',prefAngle(ivox), lineLength(ivox), lineWidth(ivox), cMap(1+floor(prefAngle(ivox)*256/pi),:))
%     hold on
% end

%order voxels by eccentricity
[eccSort indSort] = sort(roiPrf.ecc);
% [eccSort indSort] = sort(prefLevel);

% for ivox=1:length(indSort)
%     h=scatter(roiPrf.x(indSort(ivox)), roiPrf.y(indSort(ivox)),3*dotSize(indSort(ivox)),cMap(floor(prefLevel(indSort(ivox))*256/numLevels),:),'filled');
%     hold on
% end

maxSize = 0.2;
dotSize = maxSize*max(0.0001,roiPrf.r2);
dotplotsize(roiPrf.x(indSort),roiPrf.y(indSort),prefLevel(indSort),3*(dotSize(indSort).^2)./max(dotSize),30);

% h=scatter(roiPrf.x, roiPrf.y,3*dotSize,cMap(floor(prefLevel*256/numLevels),:));
% h=scatter(roiPrf.x(indSort), roiPrf.y(indSort),6*dotSize(indSort),cMap(floor(prefLevel(indSort)*256/numLevels),:));

interpImgSize = 714;
backgroundSize = 1024;
set(gca,'XTick',[]);
set(gca,'YTick',[]);
linecolor = [0.5 0.5 0.5];
line([0 0], [-backgroundSize backgroundSize],'color',linecolor);
line([-backgroundSize backgroundSize],[0 0], 'color',linecolor);
line([-interpImgSize -interpImgSize], [-interpImgSize interpImgSize], 'color',linecolor);
line([interpImgSize interpImgSize], [-interpImgSize interpImgSize], 'color',linecolor);
line([-interpImgSize interpImgSize],[interpImgSize interpImgSize], 'color',linecolor);
line([-interpImgSize interpImgSize],[-interpImgSize -interpImgSize], 'color',linecolor);
% xlim([-backgroundSize backgroundSize]); ylim([-backgroundSize backgroundSize]);
xlim([-interpImgSize interpImgSize]); ylim([-interpImgSize interpImgSize]);


axis square

end