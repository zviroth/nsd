% In brief, the nsdsynthetic experiment involved 8 runs, alternating between fixation runs and memory runs. 
% During the fixation runs, the subject performed a fixation task on a central dot; 
% during the memory runs, the subject performed a 1-back task on the presented images. 
% Each image was presented for 2 s and was followed by a 2-s gap before the next trial. 
% There are a total of 284 distinct images. 
% The stimuli and trial ordering in the 8 runs was exactly the same for all subjects (including the 1-back trials). 
% In each run, there were a total of 93 stimulus trials (as well as blank trials). 
% The 93 stimulus trials consisted of 83 regular stimulus trials plus 10 special 1-back trials (in which the presented stimulus was identical to the previously presented stimulus). 
% There were a total of 93 x 8 = 744 stimulus trials conducted in the scan session.
% There is a separate .hdf5 file created for each of the NSD subjects (based on color calibration that was tailored to each subject).
% After concatenating the 220 images with the 64 images, the result is 284 images.
% These images are shown on a gray background with RGB value (126,110,108). 
% There are 714 rows and 1360 columns. 
% The reason for the non-square shape of the image is that the word images extend far to the left and far to the right.

% <masterordering> is 1 x 744 with the sequence of trials (indices relative to the 284 images) 
% <stimpattern> is 1 session x 8 runs x 107 trials.
% elements are 0/1 indicating when stimulus trials actually occur.
close all
clear all

isub=1;

interpImgSize = 714;
backgroundSize = 1024;

%%
normResp=0;
% construct quad frequency filters
numOrientations = 4;
bandwidth = 1;
dims = [backgroundSize backgroundSize];
numLevels = maxLevel(dims,bandwidth);
[freqRespsImag, freqRespsReal, pind] = makeQuadFRs(dims, numLevels, numOrientations, bandwidth);

%%
nsdfolder = '~/NSD/';
nsdDesignFilename = fullfile(nsdfolder, 'nsd_expdesign.mat');
nsdDesign = load(nsdDesignFilename);
allImgs = nsdDesign.sharedix; %indices of the shared 1000 images
% nsdDesign.masterordering;%for each of 30000 trials, what is corresponding image (out of 10000 images)

allImgs = nsdDesign.subjectim(isub,nsdDesign.masterordering);%indices of all 10000 images used for this subject
allImgs = unique(allImgs);
%%
backgroundColor(1,1,:) = uint8([127,127,127]);

stimfolder = '~/NSD/stimuli/';
stimfilename = fullfile(stimfolder,'nsd_stimuli.hdf5');%[3 1360 714 220] 
stiminfo = h5info(stimfilename);

imgSizeX = stiminfo.Datasets.Dataspace.Size(2);%1360;
imgSizeY = stiminfo.Datasets.Dataspace.Size(3);%714;
numImgs = stiminfo.Datasets.Dataspace.Size(4);%220


nImgs = 1;

imgNum = 1;

tic
fixPoint(1,1,:) = single([255 0 0]);

% pyramidfolder = '~/NSD/stimuli/pyramid/';
pyramidfolder = ['/Volumes/nsd_sub' num2str(isub) '/pyramid/'];

iimg=0;
for imgNum=allImgs%1:1000%105:216
    iimg = iimg+1
            
    pyramidfilename = ['pyrImg' num2str(imgNum) '.mat'];
    if ~isfile(fullfile(pyramidfolder, pyramidfilename))%if file exists already no need to remake it
        origImg = h5read(stimfilename,'/imgBrick/',[1 1 1 imgNum],[3 imgSizeX imgSizeY nImgs]);
        origImg = permute(origImg,[3 2 1]);%[425,425,3]
        
        [Xq, Yq] = meshgrid(linspace(1,imgSizeX, interpImgSize), linspace(1,imgSizeY, interpImgSize));
        for irgb=1:3
            interpImg(:,:,irgb) = interp2(cast(squeeze(origImg(:,:,irgb)),'single'), Xq, Yq);
        end
        %add red semi-transparent fixation point
        
        interpImg(interpImgSize/2-8:interpImgSize/2+8,interpImgSize/2-8:interpImgSize/2+8,:) = ...
            (interpImg(interpImgSize/2-8:interpImgSize/2+8,interpImgSize/2-8:interpImgSize/2+8,:) + ...
            repmat(fixPoint,17,17,1))/2;
        
        %%add background
        bigImg = repmat(backgroundColor,backgroundSize,backgroundSize,1);
        
        bigImg(1+backgroundSize/2-interpImgSize/2 : backgroundSize/2+interpImgSize/2, 1+backgroundSize/2-interpImgSize/2 : backgroundSize/2+interpImgSize/2,:) = interpImg(:,:,:);
        
        
        % imagesc(mean(bigImg,3));
%         imagesc(bigImg);
        % imagesc(cast(bigImg,'double'));
        
        
        axis image
        % colormap gray
        
        %change to grayscale
        % for now, simply by averaging across RGB channels
        bigImg = mean(bigImg,3);
        bigImg = single(bigImg);
        
        %% pass image through steerable pyramid

        [pyr, pind] = buildQuadBands(bigImg, freqRespsImag, freqRespsReal);
        sumOri = cell(numLevels,1);
        modelOri = cell(numLevels,1);
        for ilev = 1:numLevels
            % loop over levels and orientations of the pyramid
            % initialize output
            sumOri{ilev}(:,:) = zeros(dims(1), dims(2),'single');
            modelOri{ilev} = zeros(numOrientations, dims(1), dims(2),'single');
            for orientation = 1:numOrientations
                if normResp
                    nEnergies = normEnergies(pyr,pind,numOrientations,0.1);
                    thisBand = abs(accessSteerBand(nEnergies,pind,numOrientations,ilev,orientation));
                else
                    thisBand = abs(accessSteerBand(pyr, pind, numOrientations,ilev, orientation)).^2;
                end
                sumOri{ilev}(:,:) = sumOri{ilev}(:,:) + thisBand;
                modelOri{ilev}(orientation,:,:) = thisBand;
            end
            %     sumOri{iLev}(:,:) = temp;
        end
        % for iLev = 1:numLevels
        %     levMean(iLev) = mean(abs(sumOri{iLev}(:)));
        %     levMax(iLev) = max(abs(sumOri{iLev}(:)));
        % end
        % [temp maxlev] = max(levMax);
        % maxlev
        
        save(fullfile(pyramidfolder, pyramidfilename), 'interpImgSize','backgroundSize','numOrientations',...
            'bandwidth','dims','bigImg','sumOri','modelOri','numLevels','normResp');
    end
end
toc

%% plot model output
figure
rows=1; cols=numLevels;
for ilev = 1:numLevels
    subplot(rows,cols,ilev)
    imagesc(sumOri{ilev}); 
    colormap gray
    axis image
    axis off
end
set(gcf,'position',[120 100 1200 200]);
%
figure
rows=numOrientations; cols=numLevels;
for ilev = 1:numLevels
    for orientation=1:numOrientations
    subplot(rows,cols,ilev+(orientation-1)*cols)
    imagesc(squeeze(modelOri{ilev}(orientation,:,:)));
    colormap gray
    axis image
    axis off
    end
end
set(gcf,'position',[150 50 1200 700]);