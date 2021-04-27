# LVPV
*Assessing statistical support for community detection with leading eigenvector algorithm*

This program is supplementary material to the study of the leading eigenvector algorithm for community detection \[[1](https://www.pnas.org/content/103/23/8577)\]. `lvpv.r` contains functions that allow to (i) use fine-tuning procedure in the leading algorithm; (ii) create graphs from the raw data and subsequently apply the leading eigenvector algorithm; (iii) find approximately unbiased (AU) and naive bootstrap (BP) p-values for the communities found by the leading eigenvector algorithm; and (iv) visualise results as either dendrogram or binary graph. 

The functions use numeric dataset as an input. First, a complete graph is built (via *igraph* [4] package); weights are assigned from the correlation matrix, by default an exponential transformation is used: `w = exp{(max(C)-C)/sd(C)^2}` that allows to work with negative correlation coefficients, or the coefficients can be taken into power defined by parameter beta (see *LModularity* function from give in the summary*evolQG* \[[5](https://cran.r-project.org/web/packages/evolqg/index.html)\] package for another example of such transformation for `beta = 2`). Different options for correlation are set with method and use parameters of the *cor* function. Parameter `q.cut` cuts off the smaller absolute values (with diagonal excluded) of the correlation matrix at the specified level when creating a graph. By default, a correlation matrix with any negative values is exponentially transformed. 
The leading eigenvector algorithm is applied with (`tuning = TRUE`) or without (`tuning = FALSE`) the fine-tuning procedure described in the original article \[[1](https://www.pnas.org/content/103/23/8577)\]. If parameter `bootstrap = TRUE`, Shimodaira's multiscale bootstrap procedure \[[3](https://projecteuclid.org/euclid.aos/1107794881)\] (adapted from the *pvclust* \[[2](https://cran.r-project.org/web/packages/pvclust/index.html)\] package implementation) is used to calculate approximately unbiased (AU) and naïve bootstrap (BP) estimations of the p-values for each of the detected modules. Additional features include naive bootstrap estimates for the modularity mean, standard deviation and quantiles of desired level. Results can be visualised as either graph (similar to the *pvclust* plot function) or as graph; both have an option (`members`) to highlight modules with certain significance level (or with lack of thereof) anf to label the elements. 

Additional features include naïve bootstrap estimates for the modularity mean, standard deviation and quantiles of desired level given in the summary.
See code comments for other available functions parameters.

## Examples

### Finding AU and BP p-values estimations

The bootstrap procedure and all further calculations are conducted by `PV.complete` function. It returns a list with all obtained results and used parameters. Main information can be accessed with `PV.summary` function (applied to the `PV.complete` output).
```
source("lvpv.r")

# bootstrap and p-value estimation
pvTD = PV.complete(TestData, iseed=767)
PV.summary(pvTD)

# built graph with special q.cut and diag parameters.
pvD2 = PV.complete(Data2, exp = TRUE, q.cut = .5, diag = FALSE)
PV.summary(pvD2) 

pvTD.n = PV.complete(TestData, iseed=767, exp = TRUE, q.cut = .5, diag = FALSE)
PV.summary(pvTD.n)
```

### Visualization
There are two ways to visualise results. `PV.dendro` and `PV.rectangle` functions provide dendrogram plots with annotations similar to those of the *pvclust* package, `PV.graph` with `highlight = TRUE` option highlights specified modules depending on the AU or BP values (detailes given in the code comments; see also *plot.pvclust* help). Parameter `members` allow to see the members of modules with `which` numbers (if `which` is not specified all modules are shown). 
```
# DENDROGRAM
PV.dendro(pvTD,n); PV.rectangle(pvTD.n)
# Fig. 2, right
PV.dendro(pvD2, col = "darkgreen", col.pv = c(au = "red", bp = "blue"), lwd = 2, cex.pv = 1.2, cex = 1.5)
# highlight modules with significant (e.g. AU p-value >= .95) highlighted with rectangles (Fig 2, left)
PV.dendro(pvD2); PV.rectangle(pvD2, border = "blue", fill.col = rgb(red=0, green=0, blue=1, alpha=0.1)) 
```
<img align = "center" src = "Images/dendro_options3.png" alt = "Dendrogram with rectangles" width = 850>

Information on the depicted clusters is contained within `PV.summary` output:
>Clusters:
>
>C1: a1 a2 a3 a4 a5 a6 a7  
>C2: c1 c4 c6 c7 c8  
>C3: b1 b2 b3 b4 b5  
>C4: c2 c3 c5

```
# GRAPH
PV.graph(pvTD.n) # graph with p-values
PV.graph(pvTD.n, highlight = TRUE, fill.col="lightblue") # graph with highlighted modules of AU p-value >= .95
```
In graphs, objects nodes can be hidden with `members`; labels for the objects can be fully or partially removed from the graph with `m.labels` or `which` parameters respectfully:
```
PV.graph(pvD2, highlight = TRUE, fill.col="lightblue", members = F) # hide nodes
PV.graph(pvD2, highlight = TRUE, fill.col="lightblue", m.labels = F) # hide labels (Fig. 2, left)
# show labels for the  1st and the 4th modules (Fig. 2, right):
PV.graph(pvD2, highlight = TRUE, fill.col="lightblue", which = c(1,4))
```
<img align = "center" src = "Images/graph_label_options.png" alt = "Graph labels options" width = 820>



## References
1. M.E.J. Newman, *Modularity and community structure in networks*, Proceedings of the National Academy of Sciences of the USA, 103 (2006), 8577-8582.
2. R. Suzuki, H. Shimodaira, *Pvclust: an R package for assessing the uncertainty in hierarchical
clustering*, Bioinformatics, 22 (2006), 1540–1542.
3. H. Shimodaira, *Approximately unbiased tests of regions using multistep-multiscale bootstrap
resampling*, The Annals of Statistics, 32 (2004), 2616–2641.
4. G. Csardi, T. Nepusz, *The igraph software package for complex network research*, InterJournal, Complex Systems 1695 (2006)
5. D. Melo, G. Garcia, A. Hubbe, A. P. Assis and G. Marroig, *EvolQG - An  package for evolutionary quantitative genetics*, F1000Research, 4 (2015): 925
