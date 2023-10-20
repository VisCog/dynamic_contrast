%%fMRI_joystick_analysis.m

% AmbDichop_PilotData_20230509.mat --> Contains a "results" structure with fields:
% •	subject: subject identifier
% •	hemi: hemisphere
% •	roi The Benson template predicted region of interest, where [1, 2, 3] is V1, V2, and V3, respectively. There are more values that you can find here, under “Visual Area Labels”.
% •	R2: Proportion of variance explained for each vertex.
% •	fname: bold source filename for each run, stimulus condition is in the file name
% •	bold: bold data matrix, [n_vertex, n_time], for each run

%
% NOTE: THIS IS A COPY of a script for MRI data, the actual MRI data are not saved here  


cothresh = .45;

warning('off')

subNum = 1;
ROIName = 'V1';
condNum = 1;

condList = {'LfastRslow','LslowRfast'};
condName = condList{condNum};

subName = sprintf('sub-pilot0%d',subNum);
fileName = sprintf('%s_%s_%s',subName,ROIName,condName);
load(fileName)

taskName = sprintf('task-%s',condName);
load(taskName);

%%  Choose your voxels here and average across runs

id = co>=cothresh;
tmp = zeros(size(fMRI{1}(:,id)));
for i=1:length(fMRI) % for each run
    data.experiment.response(i,:) = mean(fMRI{i}(:,id), 2);
    data.experiment.LEcontrast(i, :)  = contrastLeft;
    data.experiment.REcontrast(i, :)  =contrastRight;
end

data.t = t;

%% model fitting

freeList = {'slope', 'intercept', 'delay'};
p.delay =2; p.penalizeDelay = 4; % delay  penalization in seconds, based on previous psychophysical data
p. joystickfunction = 'delay + scale';
p.model ='b_s.softmax';
% don't let certain parameters go below zero
p.p = [1,1]; p.tau = NaN; p.m = [ 1 1];
p.k = [1,1]; p.U = [0,0,0,0];
p.sigma = 1; p.smax = 1; p.offset = 0;
freeList = {'k', 'offset'};
p.costflag = 0; p = fit('b_s.getErr', p, freeList, data);

disp(['   .. initial k left: ' num2str(round(p.k(1),3)) '   initial k right: ' num2str(round(p.k(2),3))])

% Normalize relative weights
p.k = p.k / (max(p.k));
disp(['   .. normed k left: ' num2str(round(p.k(1),3)) '    normed k right: ' num2str(round(p.k(2),3))])

% Grab model error
p.costflag = 0;  [p.step1attenuationErr,~,~,~,~,~,~] = b_s.getErr(p, data);
disp(['   .. model MSE: ' num2str(round(p.step1attenuationErr, 4)) ])
% run the fit - use k values from above, estimate U2 U3
p.U = [0 0 0 0]; freeList = {'U(2)','U(3)', 'sigma'};
%p.costflag = 1;  p = fitcon('b_s.getErr',p, freeList, dichData);
p.costflag = 1;  p = fit('b_s.getErr',p, freeList, dichData);

disp(' .. re-fitting k on all data')

p.kBeforeRefit = p.k;
freeList = {'k'};
p.costflag = 1; p = fit('b_s.getErr', p, freeList, data);
disp(['   .. initial k left: ' num2str(round(p.k(1),3)) '   initial k right: ' num2str(round(p.k(2),3))])

% Normalize relative weights
p.k = p.k / (max(p.k));
disp(['   .. normed k left: ' num2str(round(p.k(1),3)) '    normed k right: ' num2str(round(p.k(2),3))])

p.costflag = 0;
[p.softmaxErr,predModel_softmax,~,data.dich,~,~,~] = b_s.getErr(p, data);
disp(['   .. model MSE: ' num2str(round(p.softmaxErr, 4)) ])

%
%
% % Initial parameters for all model fits
% p.model = 'b_s.softmax';
% p.k = [.5,.5];
% p.p = [.25,.25];
% p.U = [0,0,0,0];
% p.sigma = .5;
% p.smax = 1;
% p.delay = 2;
% p.offset = 0;
%
% p = fitcon('getErr',p,{'delay'},stim,tSeries);
%
% freeList = {'k','U(2)>0','U(3)>0','delay','0<sigma<1','p', 'offset'};
% p = fitcon('getErr',p,freeList,stim,tSeries);
% [err,pred] = getErr(p,stim,tSeries);
% hard coded things here  :(
nReps = 3;%
preDur = 14;  %seconds
TR =1;
%% Plot model and data for whole run
figure(1)
clf
subplot('Position',[.05,.35,.9,.6])

hold on
plot(tSeries.Time,pred-min(pred),'r-');
plot(tSeries.Time,tSeries.Data-min(pred),'b-');
set(gca,'XLim',[0,max(tSeries.Time)]);
set(gca,'XLim',[min(tSeries.Time),max(tSeries.Time)]);
grid
set(gca,'XTick',[]);

legend({'pred','fMRI'});
title(sprintf('%s, err= %5.2f',fileName,err));

subplot('Position',[.05,.2,.9,.15])
image(stim.Time+p.delay,[0,1],256*stim.Data');
colormap(gray(256));
xlabel('Time (s)');
set(gca,'YTick',[0,1]);
set(gca,'YTickLabel',{'LE','RE'});
set(gca,'XLim',[min(tSeries.Time),max(tSeries.Time)]);
set(gca,'FontSize',10);

%% Plot model and data averaging across repeats within run
repLen = (max(tSeries.Time)-preDur)/(nReps*TR);  % TR's per rep

id = preDur/TR + (1:repLen);
x = tSeries.Time(id);
y = mean(reshape((tSeries.Data((preDur/TR+1):end)),repLen,nReps),2);

figure(2)
clf
subplot('Position',[.05,.35,.9,.6])
hold on
plot(x,pred(id)-min(pred),'r-');
plot(x,y-min(pred),'b-');
set(gca,'XLim',[preDur,preDur+repLen*TR]);
set(gca,'XTick',[]);
grid
legend({'pred','fMRI'});

subplot('Position',[.05,.2,.9,.15])
stimId = stim.Time>=preDur & stim.Time < preDur+repLen*TR;
image(stim.Time(stimId)+p.delay,[0,1],256*stim.Data(stimId,:)');
colormap(gray(256));
xlabel('Time (s)');
set(gca,'YTick',[0,1]);
set(gca,'YTickLabel',{'LE','RE'});
set(gca,'XLim',[preDur,preDur+repLen*TR]);
set(gca,'FontSize',10);

return
%% Plot correlation between pred and data
figure(3)
clf
plot(pred-min(pred),tSeries.Data-min(pred),'o')
axis square
xlabel('pred');
ylabel('fMRI (%)');
hold on
plot([0,max(pred)-min(pred)],[0,max(pred)-min(pred)],'r-');
grid





