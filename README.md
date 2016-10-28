
# Approach 1 (which Prof. mentioned in the 1st draft feedback) : find genes related (or statistically-significant) to two or more datasets.

**Bold faces** mean variable names / expressions in R scripts.

1.	Select **D** datasets (**D** >= 2) and name them as **all_GDS**. **D** is **length(all_GDS)**.
	- We assume in the following procedure, each dataset **all_GDS[ i ]** has:
		+ affyIDs (it's true of datasets on GPL570)
		+ **all_GDS[ i ]@dataTable@columns$disease.state** field
			- GDS5204 doesn't meet this (aging study)
			- this is for simplicity. We can generalize this later
2.	Get **common_genev** (gene ID vector common in **all_GDS**) by the following algo.
	- **affyIDs** := (all unique affy IDs in all elements of **all_GDS**)
	- translate **affyIDs** into Uniprot Accession by [this DB](https://biodbnet-abcc.ncifcrf.gov/db/db2db.php) [^1]
	- save the translation map (affy -> Uniprot) as *affy2uni.txt*
	- translate affyIDs for each dataset into Uniprot Accession by *affy2uni.txt*
	  + Note: m affyIDs <-> n Uniprot Accession
	- **common_genev** := (common genes in all datasets **all_GDS**)
	  + **length(common_genev)** will be referred to as **g** (the number of columns of the matrix M in step 3).


[^1]: Database might be incomplete. We should merge another DB if time left

# RESOURCES

+ [2-page proposal](https://docs.google.com/document/d/1WH9bjXNLgi4JiFfaLSqGhYR2SLK-xyDZ1bOGP8bEDcI/edit)

# Facts worthy of attention
+ **Some incorrect expression values in GSE are transformed into null in GDS**. Although GDS (Datasets) uses GSE (Series) as data sources, some expression values are different from GSE. In other words GDS is *"curated"* by The GEO team. So GDS is much helpful for determining which rows to be used.
+ **Some GDS needs do.log2=TRUE, and others need do.log2=FALSE.** NaN in your expression values might be due to multiple log2 transformations.
+ **getGEO() might return partial .tar.gz files without any warnings.** Using locally-cached files is strongly recommended.
+ **Some probe IDs of GPL570 and GPL96 are common.** Comparison might be possible.
