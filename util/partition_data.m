function [ytrain,ytest,ktrain] = partition_data(lY, qt, tp)
    usub = unique(qt.Submission);
    urdn = unique(qt.release_date_num);
    psub = zeros(length(urdn),1);
    for i = 1 : length(urdn)
        k = qt.release_date_num <= urdn(i);

        usubk = unique(qt.Submission(k,:));
        psub(i) = length(usubk)/length(usub);
    end

    [~,mind] = min(abs(psub - (1 - tp)));
    rdnstar = urdn(mind);
    usubtrain = unique( qt.Submission( qt.release_date_num <= rdnstar ) );
    ktrain = steq(qt.Submission, usubtrain);

    ytrain = lY(ktrain,:);
    ytest = lY(~ktrain,:);

    fprintf('DateNum cutoff: %0.0f\tDate cutoff: %s\n', rdnstar, datestr(rdnstar));
    fprintf('Samples in training set: %0.0f/%0.0f\n', size(ytrain,1), size(lY,1));
    fprintf('Samples in test set: %0.0f/%0.0f\n', size(ytest,1), size(lY,1));
    fprintf('Unique submissions in training set: %0.0f/%0.0f\n', length(usubtrain), length(usub));
    fprintf('Unique submissions in test set: %0.0f/%0.0f\n', length(usub) - length(usubtrain), length(usub));
end