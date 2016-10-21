
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

# Using Pearson's product-moment correlation coefficients in order to find genes related to two diseases.

0. Suppose we have a [ **g** x (number of samples) ] matrix *M[[ i ]]* for each dataset (1 <= i <= **D**).  
  Each entry **M[[ i ]][[ k, m ]]** is the expression value of gene **k** for sample **m**.
	- Such a matrix can be obtained from **as.data.frame(GDS2eSet(GetGEO('GDS1234')))**,  
	  whose **affyID** should be converted to **gene** later by *affy2uni.txt*.)

1. Convert the D matrices *M* into a single data frame **expdf** like:

| GDS  	| dz 	|   gene   	| eval 	|    sid    	|
|------	|:--:	|:--------:	|:----:	|:---------:	|
| 4838 	|  0 	|  X3840GE 	|  8.5 	| GSM100004 	|
| 4838 	|  0 	|  X3840GE 	|  10  	| GSM100005 	|
| 4838 	|  1 	|  X3840GE 	|   3  	| GSM100006 	|
| 4838 	|  1 	|  X3840GE 	|   1  	| GSM100007 	|
| 4838 	|  0 	| P98475AS 	|   4  	| GSM100004 	|

2. For each gene, add **cor( x = expdf$dz, y = expdf$eval  )** as a column **pearson** and the number of samples as **nsample** (including controls). The result looks like:

| GDS  	|   gene   	| pearson 	| nsample 	|
|------	|:--------:	|:-------:	|:-------:	|
| 4838 	|  X3840GE 	|   .85   	|    23   	|
| 4838 	| P98475AS 	|   .23   	|    23   	|
| 4838 	| A34578DF 	|   .42   	|    23   	|
| 2356 	|  X3840GE 	|   .88   	|    46   	|
| 2356 	| P98475AS 	|   .34   	|    46   	|

3. Select out rows with **pearson** higher than some threshold **t** (say .8). The result looks like:

| GDS  	|   gene  	| pearson 	| nsample 	|
|------	|:-------:	|:-------:	|:-------:	|
| 4838 	| X3840GE 	|   .85   	|    23   	|
| 2356 	| X3840GE 	|   .88   	|    46   	|
| 586  	| A9384AB 	|   .99   	|    4    	|

- If common genes remained for multiple rows (like X3840GE), it's worth investigating more.
- In other words: **we get interesting genes by Pearson's correlation coefficients**.  
- **CAUTION: Pearson correlation efficient for two points is always 1**. Check sample sizes and consider filtering **nsample** by some threshold *n*. Then the final result looks like:

| GDS  	|   gene  	| pearson 	| nsample 	|
|------	|:-------:	|:-------:	|:-------:	|
| 4838 	| X3840GE 	|   .85   	|    23   	|
| 2356 	| X3840GE 	|   .88   	|    46   	|


# RESOURCES

+ [2-page proposal](https://docs.google.com/document/d/1WH9bjXNLgi4JiFfaLSqGhYR2SLK-xyDZ1bOGP8bEDcI/edit)
