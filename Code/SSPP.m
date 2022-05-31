function [mapselected] = SSPP(HIM,p, mapspc,alfa,beta,sigma,~,thresold)
% Dine spatial-spectral preprocessing nonnegative matrix factorization algorithm (SSPP)
% Usage
%   [mapselected] = SSPP(HIM,p, mapspc,alfa,beta,sigma,~,thresold)
% Inputs
%   HIM - RTUD data (NC x NR x T)
%   P - the number of activity type
%   mapspc - cluster map (NC x NR)
%   alfa - percentage of spatial units with high spatial homogeneity
%   beta - percentage of spatial units with high representativeness
%   sigma - parameter for Gaussian filter
%   thresold - parameter for pixel purity index algorithm
% Outputs
%   mapselected - Abundance maps (NC x NR)
% 
% References
%   Paper: Spatial-Spectral Preprocessing Prior to Endmember Identification and Unmixing of Remotely Sensed Hyperspectral Data

[nl,ns,nb] = size(HIM);
%%%%%%% Step 1: spatial homogeneity index
h = fspecial('gaussian',15,sigma);
RMSEHIM(nl,ns) = 0;
I = imfilter(HIM,h,'symmetric');
for i = 1:nl
    for j = 1:ns
        RMSEHIM(i,j) = RMSE(HIM(i,j,:), I(i,j,:),nb);
    end
end
SpatialIndex = RMSEHIM;

%%%%%% Step 2: calcute purity for representativeness
mapSPI = SPI(HIM,p,thresold);

%%%%%%%% Step 3: select subset with high spatial homogeneity and representativeness
mapSPIselect = SelectSpectral(mapSPI,mapspc, beta);
mapSpatialselect = SelectSpatial(SpatialIndex, mapspc, alfa);
mapselected = (mapSpatialselect == mapSPIselect) .* double(mapspc);
end

function [f] = RMSE(h1,h2,nb)
f=sqrt(sum((h1-h2).^2)/nb);
end

% function for calculating purith vaule
function [mapspi] = SPI(HIM,iter,thresold)
[ns, nl, nb] = size(HIM);
[pc, ~] = pca(reshape(HIM,ns*nl,nb));
skewers = pc(:,1:iter);
mapa = zeros(ns*nl,1);
proyection = reshape(HIM,ns*nl,nb) * skewers;
for k = 1:iter   %% para k de 1 al numero de iteraciones
    maxima = max(proyection(:,k));
    minima = min(proyection(:,k));
    medio = (maxima + minima) /2;
    distmax = abs(maxima - medio);
    for i = 1:ns*nl
        peso = abs(medio - proyection(i,k)) / distmax;
        if peso > thresold
            mapa(i) = mapa(i) + peso;
        end
    end
    mapspi = reshape(mapa, ns,nl);
end
end

function [map2] = SelectSpectral(mapSPI, ClusterMap, beta)
[nl, ns] = size(mapSPI);
maxclases = max(max(ClusterMap));
map2(nl,ns) = 0;
% map2 = zeros(nl,ns,'uint8')+(maxclases+1);%%%%%%%%%
for k = 1:maxclases
    data = [0];
    count =1;
    for i = 1:nl
        for j = 1:ns
            if (ClusterMap(i,j) == k)
                data(count,1) = i;
                data(count,2) = j;
                data(count,3) = mapSPI(i,j);
                count = count +1;
            end
        end
    end
    if (count > 1)
        dataorder = sortrows(data(1:count-1,1:3),-3);
        hasta = floor((beta * count-1) / 100);
        mayoresdecero = sum(data(:,3)>0);
        if hasta > mayoresdecero
            hasta = mayoresdecero;
        end
        for l = 1:hasta
            map2(dataorder(l,1),dataorder(l,2)) = k;
        end
    end
end
end

function map2 = SelectSpatial(SpatialIndex, ClusterMap, alfa)
[nl, ns] = size(SpatialIndex);
maxclases = max(max(ClusterMap));
map2(nl,ns) = 0;
% map2 = zeros(nl,ns,'uint8')+(maxclases+2);%%%%%%%%%
for k = 1:maxclases
    data = [0];
    count =1;
    for i = 1:nl
        for j = 1:ns
            if (ClusterMap(i,j) == k)
                data(count,1) = i;
                data(count,2) = j;
                data(count,3) = SpatialIndex(i,j);
                count = count +1;
            end
        end
    end
    if (count > 1)
        dataorder = sortrows(data(1:count-1,1:3),3);
        hasta = floor((alfa * count-1) / 100);
        for l = 1:hasta
            map2(dataorder(l,1),dataorder(l,2)) = k;
        end
    end
end
end
