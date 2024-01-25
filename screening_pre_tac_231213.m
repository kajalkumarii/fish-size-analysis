clearvars,close all,clc
load('C:\PhD\data_analysis\fish-size-kinematics\231212processed-fish-size-data.mat')
%% Sanity checks
figure,hold on
histogram(fdata.gld.sp,'edgecolor','none')
histogram(fdata.brs.sp,'edgecolor','none')
xlim([0 0.25])
xlabel('Speed (m/s)')
ylabel('Number of bouts')
legend('Glide','Burst')

% figure,hold on
% plot(fdata.gld.vx,fdata.gld.vy);
% axis square

% Index only some specific stimuli:
fdata.gld.fstim = round(fdata.gld.stim-1,2);
allFlags = unique(round(fdata.gld.stim-1,2));
stims = allFlags(3:9); % Specify correct stimuli here
% stims = allFlags([3 4 5 7 8 9]); % Specify correct stimuli here

stimIdx = ismember(round(fdata.gld.stim-1,2),stims);
noCircIdx = not(fdata.gld.stim==2);
filterIdx = stimIdx & noCircIdx;

myMap = winter(1000);
blackThreshold = 15;
myMap(1:blackThreshold,:) = repmat([0 0 0],blackThreshold,1);
pts = -0.2:0.005:0.2;
% Compute density of VF:
N_VF = histcounts2(fdata.gld.vx(filterIdx), fdata.gld.vy(filterIdx), pts, pts);
N_RF = histcounts2(fdata.gld.x(filterIdx), fdata.gld.y(filterIdx), pts, pts);
figure
subplot(1,2,1)
% hist3([fdata.gld.vx(filterIdx)', fdata.gld.vy(filterIdx)'],[50 50])
imagesc(pts,pts,(N_VF)')
colormap(myMap)
cBar = colorbar
ylabel(cBar,'Number of bouts')
xlim([-.2 .2])
ylim([-.2 .2])
axis square
title('Virtual fish')
xlabel('x (m)'),ylabel('y (m)'),set(gca,'xtick',-0.2:0.1:0.2,'ytick',-0.2:0.1:0.2)
% subplot(1,2,2)
% imagesc(pts,pts,N_RF')
% set(gca,'Ydir','normal')

subplot(1,2,2)
% hist3([fdata.gld.x(filterIdx)', fdata.gld.y(filterIdx)'],[50 50])
imagesc(pts,pts,(N_RF)')
colormap(myMap)
cBar = colorbar
ylabel(cBar,'Number of bouts')
xlim([-.2 .2])
ylim([-.2 .2])
axis square
title('Real fish')
xlabel('x (m)'),ylabel('y (m)'),set(gca,'xtick',-0.2:0.1:0.2,'ytick',-0.2:0.1:0.2)

%% Check the Euclidean distance
% figure
% histogram(fdata.gld.fstim(filterIdx),'edgecolor','none')

euclDist = @(RF,VF) sqrt((RF(:,1)-VF(:,1)).^2+(RF(:,2)-VF(:,2)).^2);

distanceStims = euclDist([fdata.gld.x(filterIdx)',fdata.gld.y(filterIdx)'],...
    [fdata.gld.vx(filterIdx)',fdata.gld.vy(filterIdx)']);
figure
histogram(distanceStims,'edgecolor','none');
xlabel('Euclidean distance between VF and RF (m)')
ylabel('Number of bouts')
xlim([0 0.4])
xline(median(distanceStims),'k','linewidth',2)

distanceOverall = euclDist([fdata.gld.x',fdata.gld.y'],[fdata.gld.vx',fdata.gld.vy']);

allSizes = unique(fdata.gld.size);
allFish = unique(fdata.gld.fish);
allAges = unique(fdata.gld.age);

medStimDist = nan(numel(allSizes),numel(stims));
fishSize = nan(numel(allSizes),1);
numBouts = nan(numel(allSizes),numel(stims));
for ifish = 1:numel(allFish)
fishIdx = ismember(fdata.gld.fish,allFish(ifish)) & filterIdx;
fishSize(ifish) = unique(str2num(vertcat(fdata.gld.size{fishIdx})));
fishAge(ifish) = unique((vertcat(fdata.gld.age(fishIdx))));
    for istim = 1:numel(stims)
        stimIdx = fdata.gld.fstim==stims(istim) & fishIdx;
        medStimDist(ifish,istim) = median(distanceOverall(stimIdx));
        numBouts(ifish,istim) = sum(stimIdx);
    end
end
% Convert fishAge to number:
fishAge = fishAge';
fishAge = vertcat(fishAge{:});
fishAge = str2num(fishAge(:,2:3));

%% Sanity check: plot number of bouts per size
figure,hold on
scatter(fishSize,sum(numBouts,2),[],[0.3 0.6 0.9],'filled','markerfacealpha',0.7,'markeredgecolor','k')
xlabel('Fish size (cm)')
ylabel('Number of bouts')
model = fitlm(fishSize,sum(numBouts,2));
predictions = predict(model,[min(fishSize)-100,max(fishSize)+100]');
plot([min(fishSize)-100,max(fishSize)+100],predictions,'color',[0 0 0 0.2],'LineWidth',2)
xlim([min(fishSize)-10,max(fishSize)+10])
set(gca,'xtick',[20:10:120],'XTickLabel',[20:10:120]/100)

%% Sanity check: plot number of bouts per age
figure,hold on
scatter(fishAge,sum(numBouts,2),[],[0.3 0.6 0.9],'filled','markerfacealpha',0.7,'markeredgecolor','k')
xlabel('Fish age (dpf)')
ylabel('Number of bouts')
model = fitlm(fishAge,sum(numBouts,2));
predictions = predict(model,[min(fishAge)-100,max(fishAge)+100]');
plot([min(fishAge)-100,max(fishAge)+100],predictions,'color',[0 0 0 0.2],'LineWidth',2)
xlim([min(fishAge)-0.99,max(fishAge)+0.99])

%% Complete model
completeModel = fitlm([fishSize,fishAge],sum(numBouts,2))

%% Plotting median distances per size per fish COLORPLOT

figure
imagesc(zscore(medStimDist(:,3:end),1,2))
xlabel('Stimulus size')
ylabel('Fish size')
set(gca,"YTick",1:numel(allSizes),'YTickLabel',allSizes,"XTick",1:numel(stims),'XTickLabel',stims(3:end),'ydir','normal')
myCol = parula(1000);
colormap(flipud(myCol))
%% Plotting median distances per size per fish SCATTERPLOT (NORMALIZED)

[sortedSize,sortIdx] = sort(fishSize);
figure,hold on
sortedDist = medStimDist(sortIdx,:);
sortedCells = mat2cell(zscore(sortedDist,1,2),ones(size(sortedDist,1),1),size(sortedDist,2));
myScat = scatter(repmat(stims,1,numel(fishSize)),horzcat(sortedCells{:}),...
    [],repelem(sortedSize,numel(stims),1),'filled','markerfacealpha',0)
cBar = colorbar;
ylabel(cBar,'Fish size (cm)');
set(cBar,'ytick',[20:10:120],'yticklabel',[20:10:120]/100)
myMap = flipud(summer(numel(fishSize)));
colormap(myMap)

% Get RGB from Cdata:
Cdata = myScat.CData;
cmap = colormap;
% make it into a index image.
cmin = min(Cdata(:));
cmax = max(Cdata(:));
m = length(cmap);
index = fix((Cdata-cmin)/(cmax-cmin)*m)+1; %A
% Then to RGB
RGB = squeeze(ind2rgb(index,cmap));
RGB(:,4) = 1;
myPlot = plot(stims,zscore(sortedDist,1,2)','linewidth',2)
set(myPlot,{'Color'},num2cell(RGB(1:numel(stims):end,:),2))

xlabel('Stimulus size (cm)')
ylabel('Euclidian distance (z-score)')

%% Plotting median distances per size per fish SCATTERPLOT (ABSOLUTE)

[sortedSize,sortIdx] = sort(fishSize);
figure,hold on
sortedDist = medStimDist(sortIdx,:);
sortedCells = mat2cell(sortedDist,ones(size(sortedDist,1),1),size(sortedDist,2));
myScat = scatter(repmat(stims,1,numel(fishSize)),horzcat(sortedCells{:}),...
    [],repelem(sortedSize,numel(stims),1),'filled','markerfacealpha',0)
cBar = colorbar;
ylabel(cBar,'Fish size (cm)');
set(cBar,'ytick',[20:10:120],'yticklabel',[20:10:120]/100)
myMap = flipud(summer(numel(fishSize)));
colormap(myMap)

% Get RGB from Cdata:
Cdata = myScat.CData;
cmap = colormap;
% make it into a index image.
cmin = min(Cdata(:));
cmax = max(Cdata(:));
m = length(cmap);
index = fix((Cdata-cmin)/(cmax-cmin)*m)+1; %A
% Then to RGB
RGB = squeeze(ind2rgb(index,cmap));
RGB(:,4) = 1;
myPlot = plot(stims,sortedDist,'linewidth',2)
set(myPlot,{'Color'},num2cell(RGB(1:numel(stims):end,:),2))

xlabel('Stimulus size (cm)')
ylabel('Euclidian distance (z-score)')

%% Plotting median distances per size per fish SCATTERPLOT BINS (ABSOLUTE)

binnedDist = mat2cell(sortedDist,[6 5 5 6],7);
binnedSize = mat2cell(sortedSize,[6 5 5 6],1);
binDist = cellfun(@mean,binnedDist,'UniformOutput',false)
binSize = cellfun(@mean,binnedSize,'UniformOutput',false)
binDist = vertcat(binDist{:});
binSize = vertcat(binSize{:});

figure,hold on
myPlot = plot(stims,flipud(binDist)*100,'linewidth',2)
set(myPlot,{'Color'},num2cell(summer(4),2))
legend([num2str(round(binSize(4)/100,2)),' cm'],...
[num2str(round(binSize(3)/100,2)),' cm'],...
[num2str(round(binSize(2)/100,2)),' cm'],...
[num2str(round(binSize(1)/100,2)),' cm']);
xlabel('Stimulus size (cm)')
ylabel('Euclidian distance (cm)')


%% Plotting median distances per size per fish SCATTERPLOT BINS (NORMALIZED)

binnedDist = mat2cell(sortedDist,[6 5 5 6],7);
binnedSize = mat2cell(sortedSize,[6 5 5 6],1);
binDist = cellfun(@mean,binnedDist,'UniformOutput',false)
binSize = cellfun(@mean,binnedSize,'UniformOutput',false)
binDist = vertcat(binDist{:});
binSize = vertcat(binSize{:});

figure,hold on
myPlot = plot(stims,zscore(flipud(binDist)*100,1,2),'linewidth',2)
set(myPlot,{'Color'},num2cell(summer(4),2))
legend([num2str(round(binSize(4)/100,2)),' cm'],...
[num2str(round(binSize(3)/100,2)),' cm'],...
[num2str(round(binSize(2)/100,2)),' cm'],...
[num2str(round(binSize(1)/100,2)),' cm']);
xlabel('Stimulus size (z-score)')
ylabel('Euclidian distance (cm)')
