function d = cohens_d(group1, group2,within)
    % This function calculates the Cohen's d effect size
    % group1 - data for group/condition 1
    % group2 - data for group/condition 2
    % within - within subject? (TRUE/FALSE)     

    % Calculate the means of both groups
    mean_group1 = mean(group1);
    mean_group2 = mean(group2);

    % Calculate the standard deviations of both groups
    std_group1 = std(group1, 1);
    std_group2 = std(group2, 1);

    % Calculate the pooled standard deviation
    pooled_std = sqrt(((numel(group1) - 1) * std_group1^2 + (numel(group2) - 1) * std_group2^2) / (numel(group1) + numel(group2) - 2));

    % Compute Cohen's d
    switch within 
        case true
            dtmp = (mean_group1 - mean_group2) / pooled_std;
            r = corr(group1,group2);
            d = dtmp/sqrt(1-r);
        case false
            d = (mean_group1 - mean_group2) / pooled_std;
    end

end
