targetVar = 'Class';

selectedFeatures = { ...
    'STD', 'Quality', 'STD_QRS', ...
    'MaxIrregularityQRS3', 'MaxIrregularityQRS4', 'MaxIrregularityQRS5', ...
    'TotalIrregularityQRS', 'NumberQRS', ...
    'QRSTop3Variation', 'QRSTop3Variation2', 'QRSTop3Variation3', ...
    'pNN50', 'RMSSD', 'MedianDiffRR', 'QRSVariationX', ...
    'MeanDiffRR', 'RangeRR', 'IQR_RRraw', ...
    'SD1_SD2_ratio', ...
    'RR_Entropy', ...
    'SampleEntropy', 'LF_HF_ratio'};

minFeatures = 6;
maxFeatures = 8;
maxCombosPerK = 150;

rng(4);  % Für Reproduzierbarkeit

% --------------------------------------------------------
% FEATURE-KOMBINATIONEN SAMMELN
% --------------------------------------------------------
allCombos = {};
comboSizes = [];

for k = minFeatures:maxFeatures
    combos = nchoosek(1:length(selectedFeatures), k);

    if size(combos, 1) > maxCombosPerK
        idx = randperm(size(combos, 1), maxCombosPerK);
        combos = combos(idx, :);
    end

    for i = 1:size(combos, 1)
        allCombos{end + 1} = selectedFeatures(combos(i, :));
        comboSizes(end + 1) = k;
    end
end

numCombos = length(allCombos);
results(numCombos) = struct('Features', [], 'GMean', 0, 'Sensitivity', 0, 'Specificity', 0);

% --------------------------------------------------------
% FORTSCHRITTSANZEIGE EINRICHTEN
% --------------------------------------------------------
q = parallel.pool.DataQueue;

fprintf('Fortschritt: 0%%');
afterEach(q, @(~) updateProgress());

function updateProgress()
    persistent count
    persistent total
    if isempty(count)
        count = 0;
        total = evalin('base', 'numCombos');
    end
    count = count + 1;
    percent = (count / total) * 100;
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'); % 16 Zeichen zurück
    fprintf('Fortschritt: %3.0f%%', percent);
end

% --------------------------------------------------------
% PARALLELE AUSFÜHRUNG (SVM mit Kostenmatrix)
% --------------------------------------------------------

CostMatrix = [0 1; 15 0];  % Falsch negativ kostet 15-mal mehr

parfor i = 1:numCombos
    featSubset = allCombos{i};
    X = trainingTable(:, featSubset);
    Xmat = table2array(X);
    Y = categorical(trainingTable.(targetVar));

    try
        Mdl = fitcsvm(Xmat, Y, ...
            'KernelFunction', 'gaussian', ...
            'KernelScale', 'auto', ...
            'Standardize', true, ...
            'Cost', CostMatrix);

        cvMdl = crossval(Mdl, 'KFold', 5);
        Ypred = kfoldPredict(cvMdl);
        Ytrue = Y;

        [confMat, ~] = confusionmat(Ytrue, Ypred);

        if size(confMat, 1) == 2
            TN = confMat(1, 1);
            FP = confMat(1, 2);
            FN = confMat(2, 1);
            TP = confMat(2, 2);

            sens = TP / (TP + FN);
            spec = TN / (TN + FP);
            gmean = sqrt(sens * spec);
            if isnan(gmean), gmean = 0; end

            results(i).Features = featSubset;
            results(i).GMean = gmean;
            results(i).Sensitivity = sens;
            results(i).Specificity = spec;
        else
            results(i).Features = featSubset;
            results(i).GMean = 0;
            results(i).Sensitivity = 0;
            results(i).Specificity = 0;
        end
    catch
        results(i).Features = featSubset;
        results(i).GMean = 0;
        results(i).Sensitivity = 0;
        results(i).Specificity = 0;
    end

    send(q, i);
end

fprintf('\n');

% --------------------------------------------------------
% TOP-KOMBIS AUSGEBEN
% --------------------------------------------------------
[~, idx] = sort([results.GMean], 'descend');

topN = 10;
fprintf('\n🔝 Top %d Feature-Kombinationen nach G-Mean (√Sensitivität × Spezifität):\n', topN);
for j = 1:min(topN, numCombos)
    i = idx(j);
    fprintf('G-Mean: %.3f | Sens: %.3f | Spec: %.3f | [%s]\n', ...
        results(i).GMean, results(i).Sensitivity, results(i).Specificity, ...
        strjoin(results(i).Features, ', '));
end
