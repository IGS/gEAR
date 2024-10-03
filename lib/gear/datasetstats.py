import numpy as np
import pandas as pd

"""
Class and associated methods needed to calculate the following statistics at time of upload:
  - replicate averages
  - p-values
  - FDR (corrected p-values)
  - standard deviations
  - standard errors (of the mean)

  NOTE:
      NumPy stats: https://docs.scipy.org/doc/numpy-dev/reference/routines.statistics.html
      SciPy stats: https://docs.scipy.org/doc/scipy/reference/stats.html

      Averages: np.nanmean - https://docs.scipy.org/doc/numpy-dev/reference/generated/numpy.nanmean.html
      standard_deviations: np.nanstd - https://docs.scipy.org/doc/numpy-dev/reference/generated/numpy.nanstd.html
      standard_errors (standard error of the mean): https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.sem.html

"""


def check_all_nans(values, total_rep):
    """
    Checks to see if numpy array contains only nan values. Returns True

    Input: numpy array of replicate values

    Output: True/False
        True if all replicates are nan
        False if at least 1 value is present
    """

    all_nans = False
    if np.isnan(values).any():
        #Only do this if an nan is present.
        count_nans = 0
        for value in values:
            if np.isnan(value):
                count_nans +=1

        if count_nans == total_rep:
            all_nans = True

    return all_nans

def get_replicate_counts(obs=None):
    """
    Input: AnnData.obs observations (pandas dataframe)

    Output: python list containing the number of replicates there are for each condition.

        Example: Gene 'y' has 3 conditions. If Condition 1 has 2 replicates, Condition 2 has 3 replicates,
            and Condition 3 has 3 replicates, the list [2,3,3] is returned.
    """

    # This my 1st attempt. It worked, but it didnt utilize the power of python pandas

    # # Get replicate column
    # replicate_col = 'replicate'
    # replicates = obs.loc[:, replicate_col]
    #
    # # Build list of total number of replicates for each group
    # replicate_count = list()
    # prior_rep = None
    # for r in replicates:
    #     if prior_rep is None or r > prior_rep:
    #         prior_rep = r
    #     else:
    #         replicate_count.append(prior_rep)
    #         prior_rep = None
    #
    # # Append the last one
    # replicate_count.append(prior_rep)
    #
    # # Return the count list
    # return replicate_count

    # Get replicate counts by grouping them with pandas built-ins
    if 'time_point' in obs.columns:
        groups = obs.groupby(['cell_type', 'condition', 'time_point']).count().reset_index()
    else:
        groups = obs.groupby(['cell_type', 'condition']).count().reset_index()

    #pandas series containing number of replicates for each condition
    replicate_count = groups.loc[:, 'replicate']

    return replicate_count


class DatasetStats:
    # This class contains the methods needed to calculate various dataset statistics

    def get_replicate_averages(adata=None):
        """
        Input is AnnData object
        Output is AnnData.X np.ndarray named .Xmean

        .Xmean is the same shape as .X but contains the averaged replicate value in place of the individual replicate values
        Example:
            .X has 3 replicates for each condition (18 columns). 72, 92, 51 are 3 replicates for one condition. Their average is 71.6667.
            .Xmean also has 18 columns. 71.6667 appears 3x. Once for each replicate
            adata.X[:,0] = [ 72.  92.  51.  93.   1.  46.   0.  33.  46.  75.  56.  28.  90. 100.  7.  25.  40.  81.]
            adata.Xmean[:,0] = [71.6667 71.6667 71.6667 46.6667 46.6667 46.6667 26.3333 26.3333 26.3333 53.  53.  53.  65.6667 65.6667 65.6667 48.6667 48.6667 48.6667]
        """
        # numpy stats: https://docs.scipy.org/doc/numpy/reference/routines.statistics.html

        X = adata.X
        col_count, row_count = X.shape

        # Get replicate column
        replicate_count = get_replicate_counts(adata.obs)

        average_list = list()
        for row in range(row_count):
            replicate_vals = X[:,row]
            averages = list()
            start = 0
            end = 0
            # for total_rep in replicate_count:
            for i, total_rep in replicate_count.iteritems():
                end += total_rep

                # Subsample the replicates from 1 condition
                values = replicate_vals[start:end]

                # Is there an expression value? Skip calculation, if none found
                all_nans = check_all_nans(values, total_rep)
                if all_nans:
                    rep_avg = None
                else:
                    # Calcuate average. Ignoring NaN values and round to 4 decimals
                    rep_avg = round( np.nanmean(values), 4 )

                # Add the average for each replicate
                for r in range(total_rep):
                    averages.append(rep_avg)
                    r += 1

                # Move starting point to next set of replicates
                start += total_rep

            average_list.append(averages)
            row += 1

        # Convert list of lists to numpy array
        Xmean = np.array(average_list)
        # Xmean = pd.DataFrame(np.asarray(average_list))

        # Transpose back to match the original .X
        return Xmean.T


    def get_replicate_std(adata=None):
        """
        Input is AnnData object
        Output is AnnData.X np.ndarray named .Xstd

        .Xstd is the same shape as .X but contains the standard deviation value in place of the individual replicate values
        Example:
            .X has 3 replicates for each condition (18 columns). 72, 92, 51 are 3 replicates for one condition. Their standard deviation is 16.7398.
            .Xstd also has 18 columns. 16.7398 appears 3x. Once for each replicate
            adata.X[:,0] = [ 72.  92.  51.  93.   1.  46.   0.  33.  46.  75.  56.  28.  90. 100.  7.  25.  40.  81.]
            adata.Xstd[:,0] = [16.7398 16.7398 16.7398 37.5618 37.5618 37.5618 19.362  19.362  19.362 19.3046 19.3046 19.3046 41.684  41.684  41.684  23.669  23.669  23.669 ]
        https://docs.scipy.org/doc/numpy-dev/reference/generated/numpy.nanstd.html
        """

        X = adata.X
        col_count, row_count = X.shape

        # Get replicate column
        replicate_count = get_replicate_counts(adata.obs)

        std_list = list()
        for row in range(row_count):
            replicate_vals = X[:,row]
            stds = list()
            start = 0
            end = 0
            # for total_rep in replicate_count:
            for i, total_rep in replicate_count.iteritems():
                end += total_rep

                # Subsample the replicates from 1 condition
                values = replicate_vals[start:end]

                # Is there an expression value? Skip calculation, if none found
                all_nans = check_all_nans(values, total_rep)
                if all_nans:
                    rep_std = None
                else:
                    # Calcuate standard deviations. Ignoring NaN values and round to 4 decimals
                    rep_std = round( np.nanstd(values), 4 )

                # Add the std for each replicate
                for r in range(total_rep):
                    stds.append(rep_std)
                    r += 1

                # Move starting point to next set of replicates
                start += total_rep

            std_list.append(stds)
            row += 1

        # Convert list of lists to numpy array
        Xstd = np.array(std_list)
        # Xstd = pd.DataFrame(np.asarray(std_list))

        # Transpose back to match the original .X
        return Xstd.T


    def get_replicate_sem(adata=None):
        """
        Input is AnnData object
        Output is AnnData.X np.ndarray named .Xsem

        .Xsem is the same shape as .X but contains the standard error of the mean (SEM) in place of the individual replicate values
        Example:
            .X has 3 replicates for each condition (18 columns). 72, 92, 51 are 3 replicates for one condition. Their SEM is 11.8369.
            .Xsem also has 18 columns. 11.8369 appears 3x. Once for each replicate
            adata.X[:,0] = [ 72.  92.  51.  93.   1.  46.   0.  33.  46.  75.  56.  28.  90. 100.  7.  25.  40.  81.]
            adata.Xsem[:,0] = [11.8369 11.8369 11.8369 26.5602 26.5602 26.5602 13.691  13.691  13.691 13.6504 13.6504 13.6504 29.475  29.475  29.475  16.7365 16.7365 16.7365]
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.sem.html
        """
        import scipy.stats as stats

        X = adata.X
        col_count, row_count = X.shape

        # Get replicate column
        replicate_count = get_replicate_counts(adata.obs)

        sem_list = list()
        for row in range(row_count):
            replicate_vals = X[:,row]
            sems = list()
            start = 0
            end = 0
            # for total_rep in replicate_count:
            for i, total_rep in replicate_count.iteritems():
                end += total_rep

                # Subsample the replicates from 1 condition
                values = replicate_vals[start:end]

                # Is there an expression value? Skip calculation, if none found
                all_nans = check_all_nans(values, total_rep)
                if all_nans:
                    rep_sem = None
                else:
                    # Calcuate SEM. Ignoring NaN values and round to 4 decimals
                    rep_sem = round( stats.sem(values, nan_policy="omit"), 4 )

                # Add the SEM for each replicate
                for r in range(total_rep):
                    sems.append(rep_sem)
                    r += 1

                # Move starting point to next set of replicates
                start += total_rep

            sem_list.append(sems)
            row += 1

        # Convert list of lists to numpy array
        Xsem = np.array(sem_list)
        # Xsem = pd.DataFrame(np.asarray(sem_list))

        # Transpose back to match the original .X
        return Xsem.T


    def get_replicate_pvalue(adata=None):
        """
        Input is AnnData object
        Output is AnnData.X np.ndarray named .Xpval

        .Xpval is the same shape as .X but contains the p-value in place of the individual replicate values
        Example:
            .X has 3 replicates for each condition (18 columns).
                Control group are 3 replicates = 72, 92, 51
                Treated group are 3 replicates = 1, 46, 0
            .Xpval also has 18 columns. The pvalue to the above example 0.6836206 appears 3x. Once for each replicate.
            adata.X[:,0] = [ 72.  92.  51.  93.   1.  46.   0.  33.  46.  75.  56.  28.  90. 100.  7.  25.  40.  81.]
            adata.Xpval[:,0] = [None None None 0.6836206 0.6836206 0.6836206 None None None 0.5876048 0.5876048 0.5876048 None None None 0.4909728 0.4909728 0.4909728]

        Display for Treated when comparing Control vs treated
        One value for a set of replicates

        P-value is calculated from Pearson calculation using python scipy.stats:
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.pearsonr.html#scipy.stats.pearsonr
        """
        import scipy.stats as stats

        X = adata.X
        obs = adata.obs

        # Group conditions to get total number of replicates for each
        grouped = obs.groupby(['cell_type', 'condition', 'time_point']).count()
        grouped.reset_index().T

        # Standardized keywords that denote the control group.
        control_string = ['control', 'ctrl', 'wildtype', 'wt', 'input']

        p_values_list = list()
        for g, gene in enumerate(X.T):
            p_values = [np.nan] * len(gene)
            cell_type = None
            time_point = None
            control_list = list()
            treated_list = list()
            treated_indices = list() #track index position so p-value can later be inserted
            count_replicates = 0
            count_conditions = 1 # adds control group here

            #Go through each condition for each expression row
            for col, replicate in obs.iterrows():
                cell_type = replicate['cell_type']
                condition = replicate['condition']
                time_point = replicate['time_point']

                # Replicate is control group
                if any(condition.lower() in c for c in control_string):
                    # Add replicate to list
                    if len(control_list) < grouped.loc[cell_type, condition, time_point]['replicate']:
                        # Still adding replicates of current condition
                        control_list.append(gene[count_replicates])

                # Replicate is non-control
                else:
                    # Do we already get all the replicates for this condition?
                    if len(treated_list) < grouped.loc[cell_type, condition, time_point]['replicate']:
                        # Still adding replicates of current condition
                        treated_list.append(gene[count_replicates])
                        treated_indices.append(count_replicates)

                # All replicates accounted for. Get p-value then reset
                if len(control_list) == grouped.loc[cell_type, condition, time_point]['replicate'] and \
                    len(treated_list) == grouped.loc[cell_type, condition, time_point]['replicate']:

                    # Calculate pvalue
                    p_coef, p_value = stats.pearsonr(control_list, treated_list)
                    p_value = round(p_value, 7)

                    #Control replicates remain None, treated replicates get the pvalue
                    for t in range(len(treated_list)):
                        p_values[treated_indices[t]] = p_value

                    # Reset. Allows for additional non-control groups to be compared to the current control
                    treated_list = list()
                    treated_indices = list()
                    if count_conditions == grouped.index.levels[1].size:
                        # P-values for this set of conditions are calculated
                        # Reset control group for next condition
                        control_list = list()

                count_replicates += 1

            p_values_list.append(p_values)

        # Convert list of lists to numpy array
        Xpval = np.array(p_values_list)

        # Transpose back to match the original .X
        return pd.DataFrame(Xpval.T)


    def get_replicate_fdr(adata):
        """
        Input is AnnData object
        Output is AnnData.X np.ndarray named .Xfdr

        .Xfdr is the same shape as .X but contains the FDR (false discovery rate) in place of the individual replicate values
        Example:
            .X has 3 replicates for each condition (18 columns).
            .Xpval also has 18 columns. The pvalue to the above example 0.6836206 appears 3x. Once for each replicate.
            adata.Xpval[:,0] = [None None None 0.6836206 0.6836206 0.6836206 None None None 0.5876048 0.5876048 0.5876048 None None None 0.4909728 0.4909728 0.4909728]
            adata.Xfdr[:,0] = [None None None 0.6836206 0.6836206 0.6836206 None None None 0.6836206 0.6836206 0.6836206 None None None 0.6836206 0.6836206 0.6836206]

        FDR calculation is based on the following:
        https://stackoverflow.com/a/21739593/2900840
        """

        # Get pvalue matrix
        # Xpval = adata.Xpval.T
        Xpval = adata.uns['Xpval'].T

        # Get replicate count
        replicates = get_replicate_counts(adata.obs)

        fdrs_list = list()
        for index, row in enumerate(Xpval):
            current_fdrs = [np.nan] * len(row)
            start = 0
            end = 0
            pvalue_list = list()
            pvalue_indices = list() # tuple of int values (index_of_first_replicate, number_of_replicates)
            for i, rep_count in enumerate(replicates):
                end += rep_count

                # Take 1st replicate of group (skip None groups aka controls)
                # and store starting index and number of replicates
                if row[start:end][0] is not None:
                    pvalue_list.append(row[start:end][0])
                    pvalue_indices.append((start, rep_count))

                start += rep_count

            # The following block is heavily based on the FDR calculation provided
            #  here: https://stackoverflow.com/a/21739593/2900840
            #  orignally written by user: emre (https://stackoverflow.com/users/955278/emre)
            pvalue_list = [ (pvalue, i) for i, pvalue in enumerate(pvalue_list) ]
            pvalue_list.sort()
            pvalue_list.reverse()
            n = len(pvalue_list)
            fdr_list = list()
            for i, vals in enumerate(pvalue_list):
                rank = n - i
                pvalue, index = vals
                fdr_list.append( round((n/rank) * pvalue, 7) )
            for i in range(0, int(n)-1):
                if fdr_list[i] < fdr_list[i+1]:
                    fdr_list[i+1] = fdr_list[i]
            for i, vals in enumerate(pvalue_list):
                pvalue, index = vals
                fdr_list[index] = fdr_list[i]


            # Add the FDR to fdr_list for each replicate. Control replicates remain None
            for (i, rep_count), fdr in zip(pvalue_indices, fdr_list):
                for t in range(rep_count):
                    current_fdrs[i+t] = fdr

            fdrs_list.append(current_fdrs)

        Xfdr = np.array(fdrs_list)
        # Xfdr = pd.DataFrame(np.asarray(fdrs_list))

        # Return the transposed np.array
        return Xfdr.T
