 function [sy, nsamples] = subadjust(y, sub)
        subu = unique(sub);
        sy = zeros(size(y));
        nsamples = zeros(size(y,1),1);
        for j = 1 : length(subu)
            kk = strcmpi(sub, subu(j));
            ssize = sum(kk);
            nsamples(kk) = ssize;
            if ssize == 1
                sy(kk,:) = zeros(1, size(y,2));
            else
                %sy(kk,:) = standardize(y(kk,:));
                sy(kk,:) = bsxfun(@minus, y(kk,:), mean(y(kk,:),1));
            end
        end
    end
