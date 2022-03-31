#' LiPRawData
#'
#' Example of input LiP dataset.
#'
#' Input to MSstatsLiP converter SpectronauttoMSstatsLiPFormat.
#' Contains the following columns:
#'
#' \itemize{
#'   \item R.Condition : Label of conditions (EG Disease/Control)
#'   \item R.FileName : Name of spectral processing run
#'   \item R.Replicate : Name of biological replicate
#'   \item PG.ProteinAccessions : Protein name
#'   \item PG.ProteinGroups : Protein name, can be multiple
#'   \item PG.Quantity : Protein Quantity
#'   \item PEP.GroupingKey : Peptide grouping
#'   \item PEP.StrippedSequence : Peptide sequence
#'   \item PEP.Quantity : Peptide quantity
#'   \item EG.iRTPredicted : Predicted value
#'   \item EG.Library : Name of library
#'   \item EG.ModifiedSequence : Peptide sequence including any post-translational modifications
#'   \item EG.PrecursorId : Peptide sequence wiht modifications including charge
#'   \item EG.Qvalue : Qvalue
#'   \item FG.Charge : Identified Ion charge
#'   \item FG.Id : Peptide sequence with charge
#'   \item FG.PrecMz : Prec Mz reading
#'   \item FG.Quantity : Initial quantity reading
#'   \item F.Charge : F.Charge
#'   \item F.FrgIon : Fragment ion
#'   \item F.FrgLossType : Label for loss type
#'   \item F.FrgMz : Mz reading
#'   \item F.FrgNum : numeric Frg
#'   \item F.FrgType : character label for Frg
#'   \item F.ExcludedFromQuantification : True/False boolean for if to exclude
#'   \item F.NormalizedPeakArea : Normalized peak intensity
#'   \item F.NormalizedPeakHeight : Normalized peak height
#'   \item F.PeakArea : Unnormalized peak area
#'   \item F.PeakHeight : Unnormalized peak height
#' }
#'
#' @format A data.table consisting of 546 rows and 29 columns. Raw LiP data for use in testing
#' and examples.
#' @examples
#' head(LiPRawData)
#'
'LiPRawData'

#' TrPRawData
#'
#' Example of input TrP dataset.
#'
#' Input to MSstatsLiP converter SpectronauttoMSstatsLiPFormat.
#' Contains the following columns:
#'
#' \itemize{
#'   \item R.Condition : Label of conditions (EG Disease/Control)
#'   \item R.FileName : Name of spectral processing run
#'   \item R.Replicate : Name of biological replicate
#'   \item PG.ProteinAccessions : Protein name
#'   \item PG.ProteinGroups : Protein name, can be multiple
#'   \item PG.Quantity : Protein Quantity
#'   \item PEP.GroupingKey : Peptide grouping
#'   \item PEP.StrippedSequence : Peptide sequence
#'   \item PEP.Quantity : Peptide quantity
#'   \item EG.iRTPredicted : Predicted value
#'   \item EG.Library : Name of library
#'   \item EG.ModifiedSequence : Peptide sequence including any post-translational modifications
#'   \item EG.PrecursorId : Peptide sequence wiht modifications including charge
#'   \item EG.Qvalue : Qvalue
#'   \item FG.Charge : Identified Ion charge
#'   \item FG.Id : Peptide sequence with charge
#'   \item FG.PrecMz : Prec Mz reading
#'   \item FG.Quantity : Initial quantity reading
#'   \item F.Charge : F.Charge
#'   \item F.FrgIon : Fragment ion
#'   \item F.FrgLossType : Label for loss type
#'   \item F.FrgMz : Mz reading
#'   \item F.FrgNum : numeric Frg
#'   \item F.FrgType : character label for Frg
#'   \item F.ExcludedFromQuantification : True/False boolean for if to exclude
#'   \item F.NormalizedPeakArea : Normalized peak intensity
#'   \item F.NormalizedPeakHeight : Normalized peak height
#'   \item F.PeakArea : Unnormalized peak area
#'   \item F.PeakHeight : Unnormalized peak height
#' }
#'
#' @format A data.table consisting of 4692 rows and 29 columns. Raw TrP data for use in testing
#' and examples.
#' @examples
#' head(TrPRawData)
#'
'TrPRawData'

#' MSstatsLiP_data
#'
#' Example output of MSstatsLiP converter functions.
#'
#' Example output of MSstatsLiP converter functions. (Eg. SpectronauttoMSstatsLiPFormat).
#' A list containing two data.tables named LiP and TrP corresponding to the processed
#' LiP and TrP data now in MSstatsLiP format. The data.tables contain the following
#' columns:
#'
#' \itemize{
#'   \item ProteinName : Character column of protein names
#'   \item PeptideSequence : Character column of peptide sequence name
#'   \item PrecursorCharge : Numeric charge feature
#'   \item FragmentIon : Character fragment ion feature
#'   \item ProductCharge : Numeric charge of product
#'   \item IsotopeLabelType : Character label type
#'   \item Condition : Character label for condition (Eg. Disease/Control)
#'   \item BioReplicate : Name of biological replicate
#'   \item Run : Name of run
#'   \item Fraction : Fraction number if fractionation is present
#'   \item Intensity : Unnormalized feature intensity
#'   \item FULL_PEPTIDE(LiP data only) : Combined protein name and peptide sequence.
#'   Used for LiP data only because LiP is summarized to peptide level (not protein)
#' }
#'
#' @format A data.table consisting of 546 rows and 29 columns. Raw TrP data for use in testing
#' and examples.
#' @examples
#' head(MSstatsLiP_data$LiP)
#' head(MSstatsLiP_data$TrP)
#'
'MSstatsLiP_data'

#' MSstatsLiP_Summarized
#'
#' Example output of MSstatsLiP summarization function dataSummarizationLiP.
#'
#' Example output of MSstatsLiP summarization function dataSummarizationLiP.
#' A list containing two lists named LiP and TrP containing summarization information
#' for LiP and TrP data. Each of LiP and TrP contain data named: FeatureLevelData,
#' ProteinLevelData, SummaryMethod, ModelQC, PredictBySurvival. The two main
#' data.tables (FeatureLevelData and ProteinLevelData are shown below):
#'
#' \itemize{
#'   \item FeatureLevelData : \itemize{
#'     \item PROTEIN : Protein ID with modification site mapped in. Ex.
#'     Protein_1002_S836
#'     \item FULL_PEPTIDE (LiP Only) : Combined name of protein and peptide sequence
#'     \item PEPTIDE : Full peptide with charge
#'     \item TRANSITION: Charge
#'     \item FEATURE : Combination of Protien, Peptide, and Transition Columns
#'     \item LABEL :
#'     \item GROUP : Condition (ex. Healthy, Cancer, Time0)
#'     \item RUN : Unique ID for technical replicate of one TMT
#'     mixture.
#'     \item SUBJECT : Unique ID for biological subject.
#'     \item FRACTION : Unique Fraction ID
#'     \item originalRUN : Run name
#'     \item censored :
#'     \item INTENSITY : Original intensity value
#'     \item ABUNDANCE : Log adjusted intensity value
#'     \item newABUNDANCE : Normalized abundance column
#'   }
#'   \item ProteinLevelData : \itemize{
#'     \item RUN : MS run ID
#'     \item FULL_PEPTIDE (LiP Only) : Combined name of protein and peptide sequence
#'     \item Protein : Protein ID with modification site mapped in. Ex.
#'     Protein_1002_S836
#'     \item LogIntensities: Protein-level summarized abundance
#'     \item originalRUN : Labeling information (126, ... 131)
#'     \item GROUP : Condition (ex. Healthy, Cancer, Time0)
#'     \item SUBJECT : Unique ID for biological subject.
#'     \item TotalGroupMeasurements : Unique ID for technical replicate of one TMT
#'     mixture.
#'     \item NumMeasuredFeature : Unique ID for TMT mixture.
#'     \item MissingPercentage : Unique ID for TMT mixture.
#'     \item more50missing : Unique ID for TMT mixture.
#'     \item NumImputedFeature : Unique ID for TMT mixture.
#'   }
#'  }
#'
#' @format A list containing two lists of summarization information for LiP and TrP data.
#' @examples
#' head(MSstatsLiP_Summarized$LiP$FeatureLevelData)
#' head(MSstatsLiP_Summarized$LiP$ProteinLevelData)
#'
#' head(MSstatsLiP_Summarized$TrP$FeatureLevelData)
#' head(MSstatsLiP_Summarized$TrP$ProteinLevelData)
#'
'MSstatsLiP_Summarized'

#' MSstatsLiP_model
#'
#' Example output of groupComparisonLiP converter functions.
#'
#' Example output of MSstatsLiP groupComparisonLiP function.
#' A list containing three data.tables corresponding to unadjusted LiP, TrP, and
#' adjusted LiP models. The data.tables contain the following columns:
#'
#' \itemize{
#'   \item ProteinName : Character column of protein names
#'   \item PeptideSequence : Character column of peptide sequence name
#'   \item Label : Condition comparison (Eg. Disease vs Control)
#'   \item log2FC : Fold Change output results of model
#'   \item SE : Standard error output of model
#'   \item Tvalue : Tvalue output of model
#'   \item DF : Degrees of Freedom output of model
#'   \item pvalue : Pvalue result of model (unadjusted)
#'   \item adj.pvalue : Adjusted Pvalue, generally BH adjustement is used
#'   \item issue : Issue in model if any is reported
#'   \item MissingPercentage : Percent of missing values in specific model
#'   \item ImputationPercentage : Percent of values that needed to be imputed
#'   \item fully_TRI: Boolean indicating if Peptide is fully tryptic
#'   \item NSEMI_TRI: Boolean indicating if Peptide is NSEMI tryptic
#'   \item CSEMI_TRI: Boolean indicating if Peptide is CSEMI tryptic
#'   \item CTERMINUS: Boolean indicating if Peptide is CTERMINUS tryptic
#'   \item NTERMINUS: Boolean indicating if Peptide is NTERMINUS tryptic
#'   \item StartPos: Start position of peptide sequence
#'   \item EndPos: End position of peptide sequence
#'   \item FULL_PEPTIDE(LiP data only) : Combined protein name and peptide sequence.
#'   Used for LiP data only because LiP is summarized to peptide level (not protein)
#' }
#'
#' @format A data.table consisting of 546 rows and 29 columns. Raw TrP data for use in testing
#' and examples.
#' @examples
#' head(MSstatsLiP_model$LiP.Model)
#' head(MSstatsLiP_model$TrP.Model)
#' head(MSstatsLiP_model$Adjusted.LiP.Model)
#'
'MSstatsLiP_model'

#' SkylineTest
#'
#' Example of input data from Skylinet.
#'
#' Input to MSstatsLiP converter SkylinetoMSstatsLiPFormat
#' Contains the following columns:
#'
#' \itemize{
#'   \item Protein.Name : Name of Proteins identified by Skyline
#'   \item Peptide.Modified.Sequence : Peptide sequence
#'   \item Precursor.Charge : Charge of ion
#'   \item Fragment.Ion : Fragment ion
#'   \item Product.Charge : Identified Ion charge
#'   \item Isotope.Label.Type : Label Type
#'   \item Condition : Name of condition
#'   \item BioReplicate : name of bioreplicate annotated to data
#'   \item File.Name : Name of spectral processing run
#'   \item Area : Abudance area
#'   \item Standard.Type : Type name for row
#'   \item Truncated : Boolean if row was truncated
#' }
#'
#' @format A data.table consisting of 2115 rows and 13 columns. Raw data
#' for use in testing and examples.
#' @examples
#' head(SkylineTest)
#'
'SkylineTest'
