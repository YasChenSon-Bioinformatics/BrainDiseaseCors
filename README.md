
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
