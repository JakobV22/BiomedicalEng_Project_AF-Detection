%%
%   Viertsemesterprojekt Sommersemester S2023 
%                   - 
%   Classification in Biomedical Engineering
%
%   Author: Markus Lüken
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
    %Nach einigen Tests schein Butterwort Bandpass 2. Ordnung, 7-20Hz
    %solide zu funktionieren
    [b, a] = butter(2, [7 20]/(F_s/2), 'bandpass');
    filtered_ecg = filtfilt(b, a, ECG);
    

    

    %Wavelet Transformation
    %sym4 solide, db4 alternative
    [c,l] = wavedec(filtered_ecg, 5, 'sym4');   
    d3 = wrcoef('d', c, l, 'sym4', 3);  
    d4 = wrcoef('d', c, l, 'sym4', 4);
    d5 = wrcoef('d', c, l, 'sym4', 5);

    %Gewichtung der einzelnen Frequenzbereiche
    %4 essentiell, 3 wichtige ergänzung, 5 gute ergänzung (meistens)
    a3=1;
    a4=1;
    a5=1;
    y = a3*d3 + a4*d4 +a5*d5 ;                      %ursprünglich 4 und 5, aber 3 und 4 besser als 4 und 5
   
    %Generierung von Vergleichshöhe für Peaks
    y = abs(y).^2;
    y2 = y;
    percentileThreshold = prctile(y, 95);
    y2 = y2(y2<percentileThreshold);

    qrspeaks = [];
    
    %1. Peak-Extraction für Vergleichswert für Peak-Abstände
    [qrspeaksPRE,locsPRE] = findpeaks(y,time,'MinPeakHeight',mean(y2)*5,...
        'MinPeakDistance',0.4);
    approxPulsePRE = length(locsPRE);
    
    %2. (finale) Peak-Extraktion
    [qrspeaks,locs] = findpeaks(y,time,'MinPeakHeight',mean(y2)*5,...
        'MinPeakDistance',60/(approxPulsePRE)*0.4);      %set minDistance with knowledge about approx pulse: 
    approxPulse = length(locs);

    %find peaks that are too close and very small in local comparison
    % EFFEKTIVITÄT PRÜFEN
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
   



    %Abstände der einzelnen QRS-Komplexe
    qrsDistance = zeros(1, length(locs)-1);
    if length(locs)>1    %length(locs) mittlerweile >=38 min also egal
    for m=2:length(locs)
        qrsDistance(m-1) = locs(m)-locs(m-1);
        

    end

    else        %egal
        locs = 0;
        qrsDistance = 0;
    end

    %Standardabweichung aller QRS-Abstände
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
   
    %Paar Übliche Features
    %basierend auf einfachen Durchchnitten und Summierungen der
    %Irregularities
    noOutlierIrregularities = sortedIrregularities(4:end);    %ohne größte 5 irregularities (statt 3?)
    nNN50 = sum(noOutlierIrregularities>0.05);
    pNN50(n) = 100*nNN50/length(noOutlierIrregularities);         %pNN50 Anzahl große irregularities
    RMSSD(n) = mean(noOutlierIrregularities);                     %RMSSD durchschnittsvariation
    MedianAbsDiffRR(n) = median(qrsDistanceVariation);


    %%Rhythmusstörungen finden:

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

    %%%%%%%%%%% Schlecht zu klassifizierende class 0 fälle löschen?
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


    %weitere total scores
    QRSVariationScoreX(n) = sum(top10(4:end)); %summiert Variation scores auf
    totalIrregularitySum(n) = sum(noOutlierIrregularities); %bisschen blöd bei unterschiedlichen Pulsen
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

  %%%% QRS Distance ausreißer entfernen:
    k = 2;  % Anzahl der Extremwerte, die du entfernen willst

    qrs = qrsDistance;  % Originaldaten
    
    % Sortiere die Werte und finde die Grenzwerte
    sorted = sort(qrs);
    lowThresh = sorted(k);                  % Schwelle für kleinste k Werte
    highThresh = sorted(end - k + 1);       % Schwelle für größte k Werte
    
    % Erstelle eine Maske für alle Werte, die *nicht* zu den Extremwerten gehören
    mask = (qrs > lowThresh) & (qrs < highThresh);
    
    % Neue Version ohne die Extremwerte, Reihenfolge bleibt erhalten
    qrsFiltered = qrs(mask);

    %%%%%im folgenden qrsDistance größtenteils ersetzt durch filtered
    %%%%%variante

        % Zusätzliche Feature 1: Herzfrequenz aus RR-Intervallen
    pulseRate = 60 ./ qrsFiltered; % [bpm]
    pulseStd(n) = std(pulseRate); % Schwankung im Puls = Hinweis auf AF

    % Feature 2: mittlere Differenz zwischen RR-Intervallen
    meanDiffRR(n) = mean(abs(diff(qrsFiltered))); % absolute Änderungen zwischen Schlägen

    % Feature 3: Interquartilsabstand
    iqrRR1(n) = iqr(qrsDistance); % robustes Streuungsmaß
    iqrRR2(n) = iqr(qrsFiltered); % robustes Streuungsmaß

    % Feature 4: Maximalbereich
    rangeRR(n) = range(qrsDistance); % Max-Min der RR-Intervalle

    % Feature: Poincare -Diagramm Analyse (nichtlineare Dynamik)
    if length(qrsFiltered) > 2
    xRR = qrsFiltered(1:end-1);
    yRR = qrsFiltered(2:end);
    diffRR = yRR - xRR;
    sumRR = yRR + xRR;
    SD1(n) = std(diffRR) / sqrt(2); % kurzfristige HRV-Komponente
    SD2(n) = std(sumRR) / sqrt(2); % langfristige HRV-Komponente
    SD1_SD2_ratio(n) = SD1(n) / SD2(n); % Verhältnis als Maß für chaotische Unregelmäigkeit

    else
    SD1(n) = NaN;
    SD2(n) = NaN;
    SD1_SD2_ratio(n) = NaN;
    end

    % weitere features
    rrHist = histcounts(qrsFiltered, 20, 'Normalization', 'probability');

    rrEntropy(n) = -nansum(rrHist .* log2(rrHist + eps)); % Shannon Entropie

    rrKurtosis(n) = kurtosis(qrsFiltered); % Spitzigkeit

    rrSkewness(n) = skewness(qrsFiltered); % Schiefe



    % Einfache Approximation von Sample Entropy
    rr = qrsFiltered;
    m = 2; r = 0.2 * std(rr); N = length(rr);
    count = 0; total = 0;
    for i = 1:N - m
    for j = i+1:N - m
    if max(abs(rr(i:i+m-1) - rr(j:j+m-1))) < r
    total = total + 1;
    if abs(rr(i+m) - rr(j+m)) < r
    count = count + 1;
    end
    end
    end
    end
    if total == 0
    sampleEntropy(n) = NaN;
    else
    sampleEntropy(n) = -log(count / total + eps);
    end


    %Turning Point Ratio
    rr = qrsFiltered;
    tpr = sum((rr(2:end-1) - rr(1:end-2)) .* (rr(2:end-1) - rr(3:end)) > 0);
    tpr = tpr / (length(rr) - 2);
    turningPointRatio(n) = tpr;

    %Lempel-Ziv-Komplexität (LZC)
    rrBinary = qrsFiltered > mean(qrsFiltered);
    rrString = char(rrBinary + '0');
    uniquePatterns = {};
    i = 1;
    while i <= length(rrString)
    for j = i+1:length(rrString)
    if ~any(strcmp(rrString(i:j), uniquePatterns))
    uniquePatterns{end+1} = rrString(i:j);
    break;
    end
    end
    i = i + 1;
    end
    LZC(n) = length(uniquePatterns) / log2(length(rrBinary) + eps);



    try
    % Frequenzanalyse mit Lomb-Scargle
    rrTime = cumsum([0, qrsFiltered]); % Zeitvektor aus RR-Abständen
    [pxx, f] = plomb(qrsFiltered, rrTime(1:end-1));
    maxF = max(f); % maximale analysierte Frequenz
    % Low und High Frequency Power
    if maxF >= 0.15
    LF(n) = bandpower(pxx, f, [0.04 0.15], 'psd');
    HF(n) = bandpower(pxx, f, [0.15 0.4], 'psd');
    elseif maxF >= 0.04
    LF(n) = bandpower(pxx, f, [0.04 min(0.15, maxF)], 'psd');
    HF(n) = NaN;
    else
    LF(n) = NaN;
    HF(n) = NaN;
    end

    % Gesamte Power im Bereich
    if maxF >= 0.04
    totalPower(n) = bandpower(pxx, f, [0.04 maxF], 'psd');
    else
    totalPower(n) = NaN;
    end

    % Verhältnis LF/HF
    if ~isnan(LF(n)) && ~isnan(HF(n)) && HF(n) > 0
    LF_HF_ratio(n) = LF(n) / HF(n);
    else
    LF_HF_ratio(n) = NaN;
    end

    % Dominante Frequenz
    [~, maxIdx] = max(pxx);
    peakFreq(n) = f(maxIdx);
    catch
    % Bei Fehlern: NaN setzen
    LF(n) = NaN;
    HF(n) = NaN;
    totalPower(n) = NaN;
    LF_HF_ratio(n) = NaN;
    peakFreq(n) = NaN;
    end

%%Debug Härtefälle
if class(n) == 0 && QRSVariationScore(n) > 1
    counter = counter+1;
end
end
counter
counter
counter


%% Generate Classification Table
% Add your features to the variable names and to the table
varNames = {'STD', 'Quality', 'Class', 'STD_QRS', ...
            'MaxIrregularityQRS3', 'MaxIrregularityQRS4', 'MaxIrregularityQRS5', ...
            'TotalIrregularityQRS', 'NumberQRS', ...
            'QRSTop3Variation', 'QRSTop3Variation2', 'QRSTop3Variation3', ...
            'pNN50', 'RMSSD', 'MedianDiffRR', 'QRSVariationX', ...
            'MeanDiffRR', 'PulseStd', 'RangeRR', 'IQR_RRraw','IQR_RRfiltered', ...
            'SD1', 'SD2', 'SD1_SD2_ratio', ... 
            'RR_Entropy', 'RR_Kurtosis', 'RR_Skewness', ...
            'SampleEntropy', 'TPR', 'LZC', 'LF','HF', 'TotalPower', 'LF_HF_ratio', 'PeakFreq'};
T_class = table(standardDeviation', ...
                dataQuality', ...
                class', ...
                standardDeviationQRS', ...
                maxQRSIrregularity3', ...
                maxQRSIrregularity4', ...
                maxQRSIrregularity5', ...
                totalIrregularitySum', ...
                numberExtrQRS', ...
                QRSVariationScore', ...
                QRSVariationScore2', ...
                QRSVariationScore3', ...
                pNN50', ...
                RMSSD', ...
                MedianAbsDiffRR', ...
                QRSVariationScoreX', ...
                meanDiffRR', ...
                pulseStd', ...
                rangeRR', ...
                iqrRR1',iqrRR1', ...
                SD1', SD2', SD1_SD2_ratio', ...
                rrEntropy', rrKurtosis', rrSkewness', ...
                sampleEntropy', turningPointRatio', LZC', ...
                 LF', HF', totalPower', LF_HF_ratio', peakFreq', ...
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
      