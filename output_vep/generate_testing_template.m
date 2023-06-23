clear

responsetype = 'max';


% load random template guy
load('AM_LE_BA_15-vep-congruent.mat');

% replace fields
congruentVep.config.created = datestr(now);
congruentVep.experiment.joystickActuallyUsed = congruentVep.experiment.joystickReportOn;
congruentVep = rmfield(congruentVep, 'old_sID');
congruentVep.conditionInfo.notes{2} = 'experiment.response contains toy data for testing model fits';

switch responsetype
    case 'mean'
        congruentVep.config.sID = 'MeanExample';
        congruentVep.sID = 'MeanExample';
    case 'max'
        congruentVep.config.sID = 'MaxExample';
        congruentVep.sID = 'MaxExample';
end

% full timecourses
contrastTimecourseLE = [congruentVep.experiment.BEcontrastStart congruentVep.experiment.LEcontrast congruentVep.experiment.BEcontrastEnd];
contrastTimecourseRE = [congruentVep.experiment.BEcontrastStart congruentVep.experiment.REcontrast congruentVep.experiment.BEcontrastEnd];
contrastStacked = cat(3,contrastTimecourseLE,contrastTimecourseRE);

contrastMean = mean(contrastStacked,3);
contrastMax = max(contrastStacked,[],3);

% x-dim noise
UL = 0.5;
noiseMatrix = rescale(pinknoise(size(contrastMean)), -UL, UL);

meanPlusNoise = contrastMean + noiseMatrix;
maxPlusNoise = contrastMax + noiseMatrix(Shuffle(1:size(noiseMatrix,1)),:);

meanPlusNoise(meanPlusNoise > 1) = 1;
maxPlusNoise(maxPlusNoise > 1) = 1;

meanPlusNoise(meanPlusNoise < 0) = 0;
maxPlusNoise(maxPlusNoise < 0) = 0;

% y-dim noise
delayNoiseRange = [0.2 1.2]; % delay in sec
delayNoiseRange = round(delayNoiseRange/(1/30)); % convert to indices
delayNoise = randi(delayNoiseRange, size(congruentVep.experiment.BEcontrastStart,1),1);

responseMean = ones(size(contrastTimecourseLE))*0.5;
responseMax = responseMean;

for i = 1:size(responseMean,1)

    responseMean(i,delayNoise(i):end) =  meanPlusNoise(i,1:end-delayNoise(i)+1);
    responseMax(i,delayNoise(i):end) =  maxPlusNoise(i,1:end-delayNoise(i)+1);

end

switch responsetype
    case 'mean'
        doThisOne = responseMean;
    case 'max'
        doThisOne = responseMax;
end

tmpStr = doThisOne(:,1:size(congruentVep.experiment.BEcontrastStart,2));
tmpMid = doThisOne(:, size(congruentVep.experiment.BEcontrastStart,2)+1 : ...
    size(congruentVep.experiment.LEcontrast,2)+size(congruentVep.experiment.BEcontrastStart,2));
tmpEnd = doThisOne(:,size(congruentVep.experiment.LEcontrast,2)+size(congruentVep.experiment.BEcontrastStart,2)+1 : ...
    end);

congruentVep.experiment.binoResponseStart = tmpStr;
congruentVep.experiment.binoResponseEnd = tmpEnd;
congruentVep.experiment.response = tmpMid;

save([responsetype 'Example-vep-congruent.mat'], 'congruentVep')