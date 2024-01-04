% Joel Frohlich
% University of Tuebingen

function [lbl_out,n,statsum,idx_n,survive] = TFCE(Pvals,Tstats,dir,thresh)

% This code takes five input parameters:
%
%     Pvals: a vector of p-values Tstats: a vector of t-statistics dir: a
%     string indicating the direction of interest (positive or negative)
%     thresh: a threshold for clustering (0 --> use smith and nichols)
%     verbose: a flag indicating whether to plot the clusters (default is
%     false)
%
% The function determines whether to use a fixed or threshold-free approach
% based on the length of the thresh parameter.
%
%  The function determines which TFT elements survive based on the
%  direction of interest and the p-values and
% t-statistics. The surviving channels are then used to create a graph and
% connected components are identified.
%
% The output of the function are:
%
%     lbl_out: a vector giving each cluster a number n: the size of each
%     cluster statsum: the sum of t-statistics for each cluster idx_n: the
%     members of each cluster G: the graph of surviving elements survive:
%     the surviving elements.

% Notes If the thresh parameter is not provided, it defaults to 0.05.
% Determines the approach ('smith' or 'fixed') based on the thresh value.
% Calculates cluster mass statistics based on Mensen and Khatami (2013)
% parameters. Implements clustering based on the direction of interest,
% p-values, and t-statistics. Supports positive and negative directions.
% For threshold-free approach ('smith'), computes the support map and
% calculates cluster masses based on Smith and Nichols (2009). The function
% may plot clusters if the optional verbose parameter is set to true.

if nargin < 4
    thresh = 0.05; % clustering threshold
end

% Determine if we're using a fixed or threshold-free approach (see Smith
% and Nichols 2009, Neuroimage
% https://doi.org/10.1016/j.neuroimage.2008.03.061)
if thresh == 0
    method = 'smith';
else
    method = 'fixed';
end

% lbl_out -- gives each cluster a number n -- size of each cluster idx_n --
% members of each cluster

% parameter values below from Mensen and Khatami, 2013
E = 2/3;
H = 2;

% other values E = 1, H = 0 (cluster mass statistic) E = 0.5, H = 2 (Smith
% and Nichols, 2009) E = 1, H = 2 (Alternative approach explored by Mensen
% and Khatami, 2013)


% Determine which statistics survive thresholding

%     Initialization:
%        Create an empty matrix called survive to store information about
%        which statistical elements survive the thresholding.
% 
%     Loop through Data:
%         Go through each element of the input matrices Pvals and Tstats
%         (which contain p-values and t-statistics, respectively) using
%         nested loops.
% 
%     Direction Check:
%         Depending on the specified direction ('positive' or 'negative'),
%         check conditions for each element:
%             For 'positive' direction: Checks if the p-value is below a
%             given threshold (thresh) and if the t-statistic is greater
%             than 0. For 'negative' direction: Checks if the p-value is
%             below the threshold and if the t-statistic is less than 0.
% 
%     Update the 'survive' Matrix:
%         If the conditions are met, mark the corresponding position in the
%         survive matrix as true; otherwise, mark it as false.
% 
%     Error Handling:
%         If the specified direction is neither 'positive' nor 'negative',
%         throw an error indicating that the input is unrecognized.
% 

switch method
    case 'fixed'

        survive = nan(size(Pvals,1),size(Pvals,2));
        for i = 1:size(Pvals,1)
            for j = 1:size(Pvals,2)
                switch dir
                    case 'positive'
                        if Pvals(i,j) < thresh & Tstats(i,j) > 0 % positive direction
                            survive(ith,i,j) = true;
                        else
                            survive(ith,i,j) = false;
                        end
                    case 'negative'
                        if Pvals(i,j) < thresh & Tstats(i,j) < 0 % negative direction
                            survive(ith,i,j) = true;
                        else
                            survive(ith,i,j) = false;
                        end
                    otherwise
                        error('Input unrecognized')
                end
            end
        end

%% Cluster analysis based on threshold-free or fixed methods.

%     Surviving Elements and Connected Components:
%         We find the indices of surviving elements by converting the
%         survive matrix into logical values using logical(survive) and
%         then applying find() to get their positions. Then, we identify
%         connected components using conncomp() based on the surviving
%         elements.
% 
%     Initializing Output:
%         We set up an array called lbl_out with zeros, which will later
%         represent cluster labels.
% 
%     Checking for Empty Bins:
%         We check if there are any connected components (bins) using a
%         switch statement with isempty(bins). If there are no bins (true),
%         we set some output variables (n, idx_n, statsum) to empty arrays.
% 
%     Processing Clusters:
%         If there are bins (false), we assign cluster labels to lbl_out
%         based on the connected components. We calculate the number of
%         clusters (N) excluding the zero label, if present.
% 
%     Initializing Variables for Each Cluster:
%         We set up variables (idx_n, statsum) to store information for
%         each cluster.
% 
%     Looping Through Clusters:
%         We loop through each cluster, calculating its size (n), members
%         (idx_n), and the sum of t-statistics (statsum).
% 
%     Threshold-Free Approach (Smith's Method):
%         If the method is 'smith', we set up a range of thresholds
%         (thresh) and iterate through them. For each threshold, we
%         identify surviving elements based on the specified direction
%         ('positive' or 'negative') and the threshold.
% 
%     Updating 'survive' Matrix:
%         We mark the positions in the survive matrix as true or false
%         depending on whether the conditions are met for each element.

        survive = find(logical(survive));
        bins = conncomp(survive);

        lbl_out = zeros(1,nch);

        switch isempty(bins)
            case true
                n = [];
                idx_n = [];
                statsum = [];
            case false

                lbl_out(survive) = bins;
                if ismember(lbl_out,0)
                    N = length(unique(lbl_out))-1; % minus 1 because it also counted zero as a cluster
                else
                    N = length(unique(lbl_out));
                end
                idx_n = nan(N,nch);
                statsum = nan(N,1);

                for i = 1:N
                    n(i) = sum(lbl_out==i);
                    idx_n(i,:) = lbl_out==i;
                    statsum(i) = sum(Tstats(logical(idx_n(i,:))));
                end
        end

    case 'smith'
        thresh = linspace(0,max(abs(Tstats(:))),25);
        L = length(thresh);
        survive = nan(L,size(Pvals,1),size(Pvals,2));
        for ith = 1:L
            THR = thresh(ith);
            for i = 1:size(Pvals,1)
                for j = 1:size(Pvals,2)
                    switch dir
                        case 'positive'
                            if abs(Tstats(i,j)) > THR && Tstats(i,j) > 0 % positive direction
                                survive(ith,i,j) = true;
                            else
                                survive(ith,i,j) = false;
                            end
                        case 'negative'
                            if abs(Tstats(i,j)) > THR && Tstats(i,j) < 0 % negative direction
                                survive(ith,i,j) = true;
                            else
                                survive(ith,i,j) = false;
                            end
                        otherwise
                            error('Input unrecognized')
                    end
                end
            end
        end


%% Identify and label connected clusters of significant elements

%     Initializing the Support Map:
%         We create a 3D array named supportmap filled with NaN values,
%         which will store statistical support information for each
%         element.
% 
%     Looping Through Thresholds:
%         We iterate through a loop for each threshold value (ith) in the
%         range of thresholds (1:L).
% 
%     Extracting Surviving Elements:
%         We extract the surviving elements for the current threshold (ith)
%         from the survive matrix, transforming it into a 2D array (SRV).
% 
%     Labeling Connected Components:
%         We use bwlabel() to identify connected components (clusters) in
%         the surviving elements.
% 
%     Checking for Empty Surviving Elements:
%         We check if all elements in SRV are zero, indicating no surviving
%         elements.
% 
%     Handling Empty Case:
%         If there are no surviving elements, we initialize lbl_out with
%         zeros and set other variables (n, idx_n, and supportmap) to empty
%         values.
% 
%     Handling Non-Empty Case:
%         If there are surviving elements, we label each element in lbl_out
%         with its cluster number and determine the number of clusters (N).
% 
%     Initializing Variables:
%         We initialize variables (idx_n, stat_sum, and n) to store cluster
%         information.
% 
%     Calculating Statistics for Each Cluster:
%         We iterate through each cluster, calculating its size (n),
%         identifying its members in idx_n, and calculating the sum of
%         t-statistics (stat_sum) based on cluster members.
% 
%     Populating the Support Map:
%         We populate the supportmap with statistical support values for
%         each element based on the identified clusters.
% 
%     Summarizing Statistics:
%         We calculate statsum by squeezing and summing the supportmap
%         along the first dimension, providing an overall measure of
%         statistical support for each element in the data.

        %%%
        supportmap = nan(L,size(Pvals,1),size(Pvals,2));

        for ith = 1:L
   
            SRV = squeeze(survive(ith,:,:));
            bins = bwlabel(SRV);

            switch all(SRV(:)==0)
                case true
                    lbl_out = zeros(size(Pvals,1),size(Pvals,2));
                    n = [];
                    idx_n = [];
                    supportmap(ith,:) = 0;
                case false

                    lbl_out = bins;
                    if ismember(lbl_out,0)
                        N = length(unique(lbl_out))-1; % minus 1 because it also counted zero as a cluster
                    else
                        N = length(unique(lbl_out));
                    end

                    idx_n = nan(N,size(Pvals,1),size(Pvals,2));
                    stat_sum = nan(N,1);
                    for i = 1:N
                        n(i) = sum(lbl_out(:)==i);
                        idx_n(i,:,:) = lbl_out==i;
                        whosin = idx_n(i,:,:);
                        assert(max(whosin(:))<=1)
                        stat_sum(i) = sum(idx_n(i,:,:),[2 3])^E * thresh(ith)^H; % Eq. 1, Smith and Nichols 2009
                    end
                    for i = 1:size(Pvals,1)
                        for j = 1:size(Pvals,2)
                            if lbl_out(i,j) > 0
                                supportmap(ith,i,j) = stat_sum(lbl_out(i,j));
                            else
                                supportmap(ith,i,j) = 0;
                            end
                        end
                    end
            end

        end
        statsum = squeeze(sum(supportmap));
end

end
