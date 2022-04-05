clear;clc;
%filename = 'E:\人口分解\PopulationMixture\SampleData\RTUD_Sample.tif';
filename = '..\SampleData\RTUD_Sample.tif';
info = geotiffinfo(filename);
%sum normalization
[rtud, R] = readgeoraster(filename);
[nrow, ncol,nb] = size(rtud);
rtudnorm = rtud;
for i = 1:nrow
    for j = 1:ncol
        tmp = rtud(i,j,:);
        sumvalue = sum(tmp);
        if sumvalue >0
            rtudnorm(i,j,:) = tmp ./ sumvalue;
        end
    end
end
rtud = rtudnorm;
filecluster = '..\SampleData\RTUD_ISODATA.tif';
[cluster, R1] = readgeoraster(filecluster);
classnumber = size(unique(cluster),1);

% parameter for SSPP
actvalue = 3; %the number of activity type
alfa = 70;
beta = 30;
sigma = 1.5;
thresold = 0.7;
% SSPP method
[mapSSPP] = SSPP(rtud,actvalue,cluster,alfa,beta,sigma,classnumber,thresold);
posselect = find(mapSSPP(:));
rtudmatrix = reshape(rtud,nrow*ncol,nb);
rtud_SSPP = rtudmatrix(posselect,:);
% Temporal patterns of activities by NMF
subset = rtud_SSPP.';
colidx = find(sum(subset,1));
subset = subset(:,colidx);
option = [];
options.maxIter = 100;
options.alpha = 0;
[U_final,V_final,nIter_final,objhistory_final] = GNMF_KL(subset,actvalue,[],options); %GNMF_KL,此处执行的是标准nmf+KL散度度量方法
plot(U_final);

% activity weight paramter by constrained linear least-squares
rtudmatrix2 = rtudmatrix.';
weightMaps2 = hyperNnls(rtudmatrix2, U_final);
weightMaps2nor = weightMaps2;
for i=1:actvalue
    weightMaps2nor(i,:) = hyperNormalize(abs(weightMaps2(i,:)));
end
nSmp=size(rtudmatrix2,2);
for j=1:nSmp
    tmpsum = sum(weightMaps2nor(:,j));
    if tmpsum > 0
        weightMaps2nor(:,j) = weightMaps2nor(:,j)/(sum(weightMaps2nor(:,j)));
    end
end
%weightMap for 3d
weightMaps2nor3d = reshape(weightMaps2nor.', nrow, ncol, actvalue);

% Calculate population size for each activity
[rtud, R] = readgeoraster(filename);
rtudmatrix = reshape(rtud,nrow*ncol,nb);
number = nrow*ncol;
ratio3d = zeros(number, actvalue, nb);
PopAct3d = zeros(number, actvalue, nb);
for t = 1:nb
    for i = 1:number
        pixelw = weightMaps2nor(:,i);
        pixelw1 = pixelw.';
        pixelv = U_final(t, :);
        wv = pixelw1 .* pixelv;
        if sum(wv)>0
            ratio = wv ./sum(wv);
            ratio3d(i,:,t) = ratio;
            PopAct = ratio .* rtudmatrix(i,t);
            PopAct3d(i,:,t) = PopAct;
        end
    end
end
% output population size image
for i = 1:actvalue
    popone = PopAct3d(:,i,:);
    popone1 = reshape(popone, nrow, ncol,nb);
    route = sprintf('..\\Result\\Popsize%s.tif', num2str(i));
    geotiffwrite(route, popone1, R,'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
end

