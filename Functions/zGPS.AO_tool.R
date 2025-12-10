# ----------------------------
# Function: zGPS.AO_tool
# ----------------------------
# Purpose:
#   Perform zero–inflated Gamma–Poisson shrinker with AE ontology (zGPS.AO) analysis
#   on real drug–adverse event (AE) data to estimate group-level and AE-level relative risks.
#
# Arguments:
#   grp_data   : data.frame with columns:
#                - AE_NAME: Name of the Adverse Event (AE)
#                - GROUP_NAME: AE group name
#   pair_data  : data.frame with columns:
#                - ID: Unique identifier of the report
#                - DRUG_TYPE: Drugs reported
#                - AE_NAME: AEs reported
#                Note: If a report mentions n1 AEs and n2 drugs, generate n1*n2 rows in pair_data.
#   merge_list : Nested list specifying drugs of interest and drugs to merge
#   min_freq   : Non-negative integer for filtering scarce AEs (default 15)
#   min_size   : Positive integer for filtering small AE groups (default 20)
#
# Returns:
#   An S3 object of class 'zGPS' containing:
#   - big_data: data.frame with columns:
#       - y: Weighted count of drug-AE cell
#       - E: Baseline frequency of drug-AE cell
#       - drug: Drug type
#       - AE_grp: AE group
#       - AE: Adverse Event
#       - s: Group-level relative risk for this drug
#       - r: Overdispersion parameter for AE group
#       - p: Zero component probability for this drug-AE group
#       - beta: Log mean parameter for this drug-AE group
#       - lambda_hat: Individual AE relative risk (Bayes Estimator)
#   - pair_data: Original drug-AE pair data
#   - grp_lst: Processed AE groups after filtering
#   - merge_list: Input merge list
#   - resample_results: List of bootstrap results (initially empty)
# ----------------------------



zGPS.AO_tool = function(grp_data,
                     pair_data,
                     merge_list,
                     min_freq = 15,
                     min_size = 20) {
  grp_lst = split(grp_data, grp_data$GROUP_NAME)
  grp_lst = delete_scarse_AEs(grp_lst, pair_data, min_freq = min_freq)
  grp_lst = select_grp_size(grp_lst, min_size = min_size)
  
  big_data = suppressWarnings(get_all_params(grp_lst,
                                             pair_data,
                                             merge_list = merge_list))
  
  zGPS_result = list(
    big_data = big_data,
    pair_data = pair_data,
    grp_lst = grp_lst,
    merge_list = merge_list,
    resample_results = list()
  )
  class(zGPS_result) = 'zGPS'
  
  return(zGPS_result)
}
