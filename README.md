# MetFlowChord
We use the circlize[1] package to draw the following chord graphs with [R scripts](https://github.com/biowilliam/MetFlowChord/blob/master/Chord%20visualization%20Code%202020.R)

1. Untreated PBMC celltype markers vs Metabolic markes plot: 
The spearman correlation value was calculated by the gMFI between each metabolic protein and each immune cell subset markers.
It is between a cell type vs all others ( e.g. T CD8+ vs All Others cells ) by gMFI of each marker. Negative values of gMFI are replaced as 0 in the dataset  (gMFI_data_Untreated.xlsx) before the spearman calculation. The

 ![alt text](https://github.com/biowilliam/MetFlowChord/blob/master/Untreated%20PBMC.png)
 
2. T cell Immuno vs Metabolic markers plot: 
Calculate the difference of CD4+ gMFI between activated and untreated T cells for each metabolic protein and each immune marker. Then we get the the spearman correlation between those metabolic and immune markers.

![alt text](https://github.com/biowilliam/MetFlowChord/blob/master/T%20cell%20Immuno%20vs%20Metabolic.png)


## Sources
 1. [Zuguang Gu, Lei Gu, Roland Eils, Matthias Schlesner, Benedikt Brors, circlize Implements and enhances circular visualization in R. Bioinformatics (Oxford, England) 2014.](https://academic.oup.com/bioinformatics/article/30/19/2811/2422259)
