%%
%   Viertsemesterprojekt Sommersemester S2023 
%                   - 
%   Classification in Biomedical Engineering
%
%   Author: Markus L�ken
%   Email: lueken@hia.rwth-aachen.de
%
%%
close all;
clear all;


F_s = 200; % [Hz]
scale = 0.0011; % [mV/bit]

% If you like to plot the data, set to one
plotECGdata = 0;

%% Load the data overview file
% Keep this code unchanged
datTable = readtable('ECGs_CSV_train\_DataFile.csv');


%%%%%%%%%%
load("TrainingTable.mat");
trQRSVariationScore1 = trainingTable.QRSTop3Variation;
sortedTrQRS1 = sort(trQRSVariationScore1, 'descend');
sortedTrQRS1(50)


counter=0;
%%%%%%%%%%
%% Process Data Files (Feature Extraction)
for n=1:height(datTable)
    disp(['Currently processing file ', datTable.fileName{n},...
          ', ', num2str(n/height(datTable)*100) '% (', num2str(n),...
          '/' , num2str(height(datTable)), ') completed.']);
   
    data_raw = load(['ECGs_CSV_train\' datTable.fileName{n}]); 
   
    time = ([1:length(data_raw)]-1)/F_s;
    ECG = data_raw*scale;

    
    
    class(n) = datTable.doc_eval_0SR_1AF_2Other_3Uninterpretable(n);
    if class(n) == 3
        class(n) = [];
        continue;
    end
    %%
    %
    % Your code goes here


    locs = [];
   
    
    %Filter?
    [b, a] = butter(2, [7 20]/(F_s/2), 'bandpass');
    filtered_ecg = filtfilt(b, a, ECG);
    

    

    %TEST
    %anderer wavelet ansatz chatgpt
    [c,l] = wavedec(filtered_ecg, 5, 'sym4');   %db4 bisschen besser als sym4 aktuell
    d3 = wrcoef('d', c, l, 'sym4', 3);  
    d4 = wrcoef('d', c, l, 'sym4', 4);
    d5 = wrcoef('d', c, l, 'sym4', 5);
    a3=1;
    a4=1;
    a5=1;
    y = a3*d3 + a4*d4 +a5*d5 ;                      %urspr�nglich 4 und 5, aber 3 und 4 besser als 4 und 5
   
    y = abs(y).^2;
    y2 = y;
    percentileThreshold = prctile(y, 95);
    y2 = y2(y2<percentileThreshold);
    qrspeaks = [];
    
    [qrspeaksPRE,locsPRE] = findpeaks(y,time,'MinPeakHeight',mean(y2)*5,...
        'MinPeakDistance',0.4);
    approxPulsePRE = length(locsPRE);
    %%%TEST
    %%approxPeakHeight = median(qrspeaksPRE);
    %%%TEST
    [qrspeaks,locs] = findpeaks(y,time,'MinPeakHeight',mean(y2)*5,...
        'MinPeakDistance',60/(approxPulsePRE)*0.4);      %set minDistance with knowledge about approx pulse
    approxPulse = length(locs);

    %find peaks that are too close and very small in local comparison
    % EFFEKTIVIT�T PR�FEN
    %irgendwas klappt nicht
    k=2;
    while k< length(locs)-1

        locDif1 = locs(k)-locs(k-1);
        locDif2 = locs(k+1)-locs(k);
        locDifRelation = locDif1/locDif2;

        if  (locDif1 + locDif2 < (120/approxPulse) *0.75)||(locDif1< (60/approxPulse*0.5))||(locDif2< (60/approxPulse*0.5))||(locDifRelation >1.5) || (locDifRelation <0.66)
            if ((qrspeaks(k)<(qrspeaks(k+1)/1.5)) && (qrspeaks(k)<(qrspeaks(k-1)/1.5)))

                if ((qrspeaks(k)<(qrspeaks(k+1)/3)) || (qrspeaks(k)<(qrspeaks(k-1)/3)))
                qrspeaks(k) = [];
                locs(k) = [];
                continue;

                end
            end
        end
        k=k+1;
    end
    qrspeaks(1) = [];
    locs(1) = [];
    qrspeaks(end) = [];
    locs(end) = [];
   



    %Abst�nde der einzelnen QRS-Komplexe
    qrsDistance = zeros(1, length(locs)-1);
    if length(locs)>1    %locs mittlerweile >=38 min also egal
    for m=2:length(locs)
        qrsDistance(m-1) = locs(m)-locs(m-1);
        

    end

    else        %egal
        locs = 0;
        qrsDistance = 0;
    end

    %Standardabweichung aller QRS-Abst�nde
    standardDeviationQRS(n) = std(qrsDistance);
    
    %Find Differences between adjacent qrs intervals
     
    if length(qrsDistance)>1
    qrsDistanceVariation = zeros(1, length(qrsDistance)-1);
    for m=2:length(qrsDistance)
    qrsDistanceVariation(m-1) = abs(qrsDistance(m)-qrsDistance(m-1));
    end
    else
        qrsDistanceVariation = 0;
    end

    %Find max irregularities (many likely caused by misdetection)
    sortedIrregularities = sort(qrsDistanceVariation, 'descend');
    maxQRSIrregularity(n) = sortedIrregularities(1);
    maxQRSIrregularity2(n) = sortedIrregularities(2);
    maxQRSIrregularity3(n) = sortedIrregularities(3);
    maxQRSIrregularity4(n) = sortedIrregularities(4);
    maxQRSIrregularity5(n) = sortedIrregularities(5);
   
    %Paar �bliche Features
    noOutlierIrregularities = sortedIrregularities(4:end);    %ohne gr��te 5 irregularities (statt 3?)
    nNN50 = sum(noOutlierIrregularities>0.05);
    pNN50(n) = 100*nNN50/length(noOutlierIrregularities);         %pNN50 Anzahl gro�e irregularities
    RMSSD(n) = mean(noOutlierIrregularities);                     %RMSSD durchschnittsvariation
    MedianAbsDiffRR(n) = median(qrsDistanceVariation);


    %%Rhythmusst�rungen finden:

    top10 = [0,0,0,0,0,0,0,0,0,0];
    top10locs = [0,0,0,0,0,0,0,0,0,0];
    for m=1:1:length(qrsDistance)-4
        window = [qrsDistance(m),qrsDistance(m+1),qrsDistance(m+2),qrsDistance(m+3)];
        windowMedian = median(window);
        
        windowVariations = abs(window -windowMedian);
        variationScore = sum(windowVariations);
        for l=1:10
            if variationScore>top10(l)
                top10(l) = variationScore;
                top10locs(l) = m;
                break;
            end
        end
    end

    QRSVariationScore(n) = top10(4);
    QRSVariationScore2(n) = top10(6);
    QRSVariationScore3(n) = top10(8);

    %%%%%%%%%%% Schlecht zu klassifizierende class 0 f�lle l�schen?
    % if class(n) == 0 && QRSVariationScore3(n) > 1
    %     QRSVariationScore(n) = [];
    %     QRSVariationScore2(n) = [];
    %     QRSVariationScore3(n) = [];
    %     pNN50(n) =[];
    %     RMSSD(n) = [];
    %     MedianAbsDiffRR(n) = [];
    %     maxQRSIrregularity(n) = [];
    %     maxQRSIrregularity2(n) = [];
    %     maxQRSIrregularity3(n) = [];
    %     maxQRSIrregularity4(n) = [];
    %     maxQRSIrregularity5(n) = [];
    %     standardDeviationQRS(n) = [];
    %     class(n) = [];
    % 
    %     continue;
    % end
    %%%%%%%%%%
    QRSVariationScoreX(n) = sum(top10(4:end));
    totalIrregularitySum(n) = sum(qrsDistanceVariation);
    length(locs)
    numberExtrQRS(n) = length(locs);
    %
    % Example:
    standardDeviation(n) = std(ECG);
    dataQuality(n) = datTable.signalQuality_0Excellent_1Good_2Bad_3Uninterpretable(n);
    %
    %%

    % Plot data for visual inspection
    if plotECGdata==1 && class(n) == 0 && QRSVariationScore3(n) > 1
    f1 = figure('Visible','on','Name','ECG Plot','NumberTitle','off');
    plot(time,y)
    hold on
    plot(locs,qrspeaks,'ro')
    xlabel('Seconds')
    title('R Peaks Localized by Wavelet Transform with Automatic Annotations, Pulse:', length(locs))
   
    %top10locs(1)
    %top10locs(2)
    %top10locs(3)
               %
    f2 = figure('Name', 'ECG Raw', 'NumberTitle', 'off');
    plot(time, ECG)
   
    ylabel('ECG [mV]')
    xlabel('Time [s]')
    if datTable.doc_eval_0SR_1AF_2Other_3Uninterpretable(n)==0
        title([datTable.fileName{n}, ': Sinus Rhythm' ])
    elseif datTable.doc_eval_0SR_1AF_2Other_3Uninterpretable(n)==1
        title([datTable.fileName{n}, ': Atrial Fibrillation' ])
    elseif datTable.doc_eval_0SR_1AF_2Other_3Uninterpretable(n)==3
        title([datTable.fileName{n}, ': Uninterpretable' ])
    end
    
    f3 = figure('Name', 'ECG Filtered', 'NumberTitle', 'off');
    plot(time, filtered_ecg)
   
    ylabel('ECG [mV]')
    xlabel('Time [s]')
        
    %pause;
    waitfor(f1)
    waitfor(f2)
    waitfor(f3)
    end
if class(n) == 0 && QRSVariationScore(n) > 1
    counter = counter+1;
end
end
counter
counter
counter


%% Generate Classification Table
% Add your features to the variable names and to the table
varNames = {'STD',...
            'Quality',...
            'Class', 'STD_QRS', 'MaxIrregularityQRS3', 'MaxIrregularityQRS4', 'MaxIrregularityQRS5' 'TotalIrregularityQRS', 'NumberQRS', 'QRSTop3Variation', 'QRSTop3Variation2', 'QRSTop3Variation3', 'pNN50', 'RMSSD', 'MedianDiffRR',...
            'QRSVariationX'};

T_class = table(standardDeviation',...
                dataQuality',class',standardDeviationQRS',maxQRSIrregularity3',maxQRSIrregularity4',maxQRSIrregularity5',...
                totalIrregularitySum', numberExtrQRS',QRSVariationScore',QRSVariationScore2',QRSVariationScore3',pNN50', RMSSD', MedianAbsDiffRR',QRSVariationScoreX',...
                'VariableNames', varNames);

            
% save the training data
trainingTable = T_class;

save('TrainingTable.mat', 'trainingTable');
      
      
%% Run the Classification Learner and save your final classification model via 'Export -> Generate Function' --> save as 'trainClassifier'    

%% Test the classification learning and accuracy
load('TrainingTable.mat')

[myClassifier, validationAccuracy] = trainClassifier(trainingTable);
disp(['Classification Accuracy after Crossvalidation: ' num2str(validationAccuracy, '%1.3f')]);

% Classify the training data (later, the test data will be classified here)
yfit = myClassifier.predictFcn(trainingTable);    
      