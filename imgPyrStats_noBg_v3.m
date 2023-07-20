close all
clear all

reloadData = 1;

nsdFolder = '~/NSD/stimuli/';
nsdFolder = '/misc/data18/rothzn/nsd/stimuli/'; %save file
% pyramidfolder = '/System/Volumes/Data/misc/data18/rothzn/nsd/stimuli/pyramid';
pyramidfolder = '/misc/data18/rothzn/nsd/stimuli/pyramid_noBackground_v3';
interpImgSize = 714;
backgroundSize = 1024;
imgScaling = 0.5;
numOrientations = 8;
bandwidth = 1;
dims = [backgroundSize backgroundSize];
dims = dims*imgScaling;
numLevels = maxLevel(dims,bandwidth);

nimgs = 73000;
% nimgs = 1000;
if reloadData
    tic
    for ilev=1:numLevels
        allModelOri{ilev} = zeros(numOrientations,dims(1),dims(2));
    end
    for imgNum = 1:nimgs
        imgNum
        pyramidfilename = ['pyrImg' num2str(imgNum) '.mat'];
        load(fullfile(pyramidfolder, pyramidfilename), 'interpImgSize','backgroundSize','imgScaling',...
            'numOrientations','bandwidth','dims','bigImg','sumOri','modelOri','numLevels','normResp');

        for ilev = 1:numLevels
            for orientation = 1:numOrientations
                %average model output across images for each SF and orientation filter
                allModelOri{ilev}(orientation,:,:) = allModelOri{ilev}(orientation,:,:) + modelOri{ilev}(orientation,:,:)./nimgs;

            end
        end
    end

    save([nsdFolder 'imgPyrStats_noBg_v3_' num2str(nimgs) 'imgs.mat'],"allModelOri",'nimgs','numLevels','numOrientations');
    toc
else
    load([nsdFolder 'imgPyrStats_noBg_v3_' num2str(nimgs) 'imgs.mat'],"allModelOri",'nimgs','numLevels','numOrientations');
end

figure(1);
rows=numLevels;
cols = numOrientations;
isubplot=0;
for ilev = 1:numLevels
    for orientation = 1:numOrientations
        isubplot=isubplot+1;
        subplot(rows,cols,isubplot)
        imagesc(squeeze(allModelOri{ilev}(orientation,:,:)));
        axis off
        axis square
    end
end
set(gcf,'position',[400 100 1100 900])