```R

```


```R
# Data used
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79761
# download raw data
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE79761&format=file
```


```R

```


```R
options(warn=-1)
```


```R
#load required libraries
library(affy)
library(limma)
library(oligo)
```


```R
setwd(R"(D:\tutorial\GSE79761_RAW)") #set the working directory
```


```R
#read files with *.cel file patterns
affyRaw <- read.celfiles(list.files(pattern = '*CEL.gz', full.names = TRUE))
```

    Platform design info loaded.
    
    

    Reading in : ./GSM2102086_1V.CEL.gz
    Reading in : ./GSM2102087_2V.CEL.gz
    Reading in : ./GSM2102088_3V.CEL.gz
    Reading in : ./GSM2102089_1E2.CEL.gz
    Reading in : ./GSM2102090_2E2.CEL.gz
    Reading in : ./GSM2102091_3E2.CEL.gz
    Reading in : ./GSM2102092_1D.CEL.gz
    Reading in : ./GSM2102093_2D.CEL.gz
    Reading in : ./GSM2102094_3D.CEL.gz
    Reading in : ./GSM2102095_1DE.CEL.gz
    Reading in : ./GSM2102096_2DE.CEL.gz
    Reading in : ./GSM2102097_3DE.CEL.gz
    


```R
#normaliza data
eset = rma(affyRaw, normalize=TRUE)
```

    Background correcting
    Normalizing
    Calculating Expression
    


```R
#write the normalized value to the local system
write.exprs(eset, file = "GSE79761_rma_normalized.txt")
```


```R
rownames(eset)[1:10] #view 10 rownames
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'1007_s_at'</li><li>'1053_at'</li><li>'117_at'</li><li>'121_at'</li><li>'1255_g_at'</li><li>'1294_at'</li><li>'1316_at'</li><li>'1320_at'</li><li>'1405_i_at'</li><li>'1431_at'</li></ol>




```R
# library(biomaRt)
```


```R
# mart <- useMart("ENSEMBL_MART_ENSEMBL")
# mart <- useDataset("hsapiens_gene_ensembl", mart)
```


```R
# annotLookup <- getBM(mart = mart,
#                      attributes = c("affy_hg_u133_plus_2", "ensembl_gene_id", "gene_biotype", "external_gene_name"),
#                      filter = "affy_hg_u133_plus_2", values = rownames(eset), uniqueRows=TRUE)
```


```R
# head(annotLookup)
```


```R
#load the targets: sample information
targets = read.table("targets.txt", header=T, sep='\t', row.names = 1)
head(targets,12)
```


<table class="dataframe">
<caption>A data.frame: 12 × 3</caption>
<thead>
	<tr><th></th><th scope=col>CellLine</th><th scope=col>Group</th><th scope=col>Replicate</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>GSM2102086</th><td>MCF7</td><td>Veh   </td><td>Rep1</td></tr>
	<tr><th scope=row>GSM2102087</th><td>MCF7</td><td>Veh   </td><td>Rep2</td></tr>
	<tr><th scope=row>GSM2102088</th><td>MCF7</td><td>Veh   </td><td>Rep3</td></tr>
	<tr><th scope=row>GSM2102089</th><td>MCF7</td><td>Est   </td><td>Rep1</td></tr>
	<tr><th scope=row>GSM2102090</th><td>MCF7</td><td>Est   </td><td>Rep2</td></tr>
	<tr><th scope=row>GSM2102091</th><td>MCF7</td><td>Est   </td><td>Rep3</td></tr>
	<tr><th scope=row>GSM2102092</th><td>MCF7</td><td>Dex   </td><td>Rep1</td></tr>
	<tr><th scope=row>GSM2102093</th><td>MCF7</td><td>Dex   </td><td>Rep2</td></tr>
	<tr><th scope=row>GSM2102094</th><td>MCF7</td><td>Dex   </td><td>Rep3</td></tr>
	<tr><th scope=row>GSM2102095</th><td>MCF7</td><td>DexEst</td><td>Rep1</td></tr>
	<tr><th scope=row>GSM2102096</th><td>MCF7</td><td>DexEst</td><td>Rep2</td></tr>
	<tr><th scope=row>GSM2102097</th><td>MCF7</td><td>DexEst</td><td>Rep3</td></tr>
</tbody>
</table>




```R
# eset@assayData$exprs
```


```R
#read the normalized data saved in local 
df = read.table('GSE79761_rma_normalized.txt', sep='\t', header = T, row.names = 1)
dim(df)
head(df,4)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>54675</li><li>12</li></ol>




<table class="dataframe">
<caption>A data.frame: 4 × 12</caption>
<thead>
	<tr><th></th><th scope=col>GSM2102086_1V.CEL.gz</th><th scope=col>GSM2102087_2V.CEL.gz</th><th scope=col>GSM2102088_3V.CEL.gz</th><th scope=col>GSM2102089_1E2.CEL.gz</th><th scope=col>GSM2102090_2E2.CEL.gz</th><th scope=col>GSM2102091_3E2.CEL.gz</th><th scope=col>GSM2102092_1D.CEL.gz</th><th scope=col>GSM2102093_2D.CEL.gz</th><th scope=col>GSM2102094_3D.CEL.gz</th><th scope=col>GSM2102095_1DE.CEL.gz</th><th scope=col>GSM2102096_2DE.CEL.gz</th><th scope=col>GSM2102097_3DE.CEL.gz</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1007_s_at</th><td>11.714877</td><td>12.485145</td><td>11.332145</td><td>11.629958</td><td>11.358790</td><td>11.287053</td><td>11.660072</td><td>11.379074</td><td>11.391383</td><td>11.570162</td><td>12.359483</td><td>11.482696</td></tr>
	<tr><th scope=row>1053_at</th><td> 8.985365</td><td> 9.283224</td><td> 9.092368</td><td> 9.392458</td><td> 9.198721</td><td> 9.236350</td><td> 9.262318</td><td> 9.041983</td><td> 8.994891</td><td> 9.505972</td><td> 9.265597</td><td> 9.078676</td></tr>
	<tr><th scope=row>117_at</th><td> 6.570190</td><td> 6.053559</td><td> 6.318432</td><td> 6.569726</td><td> 6.400025</td><td> 6.400079</td><td> 6.428471</td><td> 5.945971</td><td> 6.100483</td><td> 6.449865</td><td> 5.946788</td><td> 6.378356</td></tr>
	<tr><th scope=row>121_at</th><td> 8.230666</td><td> 8.370094</td><td> 8.178881</td><td> 8.255839</td><td> 8.172376</td><td> 8.178926</td><td> 8.260911</td><td> 8.143419</td><td> 8.218729</td><td> 8.209071</td><td> 8.251734</td><td> 8.218064</td></tr>
</tbody>
</table>




```R
#create design matrix based on the group
conditions<- paste(targets$Group)
head(conditions)
conditions <- factor(conditions, levels=unique(conditions))
#create design matrix
design <- model.matrix(~0+ conditions)
colnames(design) <- levels(conditions) #add column name to design matrix
unique(conditions)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Veh'</li><li>'Veh'</li><li>'Veh'</li><li>'Est'</li><li>'Est'</li><li>'Est'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>Veh</li><li>Est</li><li>Dex</li><li>DexEst</li></ol>

<details>
	<summary style=display:list-item;cursor:pointer>
		<strong>Levels</strong>:
	</summary>
	<style>
	.list-inline {list-style: none; margin:0; padding: 0}
	.list-inline>li {display: inline-block}
	.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
	</style>
	<ol class=list-inline><li>'Veh'</li><li>'Est'</li><li>'Dex'</li><li>'DexEst'</li></ol>
</details>



```R
#view design matrix
design
```


<table class="dataframe">
<caption>A matrix: 12 × 4 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>Veh</th><th scope=col>Est</th><th scope=col>Dex</th><th scope=col>DexEst</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>1</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>2</th><td>1</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>3</th><td>1</td><td>0</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>4</th><td>0</td><td>1</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>5</th><td>0</td><td>1</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>6</th><td>0</td><td>1</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>7</th><td>0</td><td>0</td><td>1</td><td>0</td></tr>
	<tr><th scope=row>8</th><td>0</td><td>0</td><td>1</td><td>0</td></tr>
	<tr><th scope=row>9</th><td>0</td><td>0</td><td>1</td><td>0</td></tr>
	<tr><th scope=row>10</th><td>0</td><td>0</td><td>0</td><td>1</td></tr>
	<tr><th scope=row>11</th><td>0</td><td>0</td><td>0</td><td>1</td></tr>
	<tr><th scope=row>12</th><td>0</td><td>0</td><td>0</td><td>1</td></tr>
</tbody>
</table>




```R
#fit the design to the expression data
fit = lmFit(df, design)
#make contrast matrix
contrast.matrix <- makeContrasts(degs = Dex - Veh,
                                 degs2 = Est - Veh,
                                 levels = design)
#Fit contrast matrix
fit.cont = contrasts.fit(fit, contrast.matrix)
#Fit ebayes
fit.cont = eBayes(fit.cont)
```


```R
#view contrast matrix
contrast.matrix
```


<table class="dataframe">
<caption>A matrix: 4 × 2 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>degs</th><th scope=col>degs2</th></tr>
</thead>
<tbody>
	<tr><th scope=row>Veh</th><td>-1</td><td>-1</td></tr>
	<tr><th scope=row>Est</th><td> 0</td><td> 1</td></tr>
	<tr><th scope=row>Dex</th><td> 1</td><td> 0</td></tr>
	<tr><th scope=row>DexEst</th><td> 0</td><td> 0</td></tr>
</tbody>
</table>




```R
#generate the degs 
topTable(fit.cont, coef = "degs", adjust.method = "BH", number = nrow(df), sort.by = 'p')
```


<table class="dataframe">
<caption>A data.frame: 54675 × 6</caption>
<thead>
	<tr><th></th><th scope=col>logFC</th><th scope=col>AveExpr</th><th scope=col>t</th><th scope=col>P.Value</th><th scope=col>adj.P.Val</th><th scope=col>B</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>212448_at</th><td> 0.7952617</td><td> 7.715136</td><td> 7.515949</td><td>1.692823e-05</td><td>0.8484434</td><td>-2.835903</td></tr>
	<tr><th scope=row>228865_at</th><td> 1.0253649</td><td> 5.512106</td><td> 6.288833</td><td>7.860038e-05</td><td>0.8484434</td><td>-2.980562</td></tr>
	<tr><th scope=row>203973_s_at</th><td> 1.1970560</td><td>10.202935</td><td> 6.015100</td><td>1.136349e-04</td><td>0.8484434</td><td>-3.021213</td></tr>
	<tr><th scope=row>241782_at</th><td> 1.2243888</td><td> 4.757880</td><td> 5.930289</td><td>1.276446e-04</td><td>0.8484434</td><td>-3.034554</td></tr>
	<tr><th scope=row>203542_s_at</th><td> 1.9430377</td><td> 9.142446</td><td> 5.924000</td><td>1.287549e-04</td><td>0.8484434</td><td>-3.035558</td></tr>
	<tr><th scope=row>227703_s_at</th><td> 1.3467424</td><td> 8.883694</td><td> 5.898666</td><td>1.333332e-04</td><td>0.8484434</td><td>-3.039624</td></tr>
	<tr><th scope=row>230636_s_at</th><td> 1.2840063</td><td> 6.532282</td><td> 5.797603</td><td>1.534103e-04</td><td>0.8484434</td><td>-3.056179</td></tr>
	<tr><th scope=row>223168_at</th><td> 1.6276482</td><td> 9.938139</td><td> 5.638520</td><td>1.918613e-04</td><td>0.8484434</td><td>-3.083369</td></tr>
	<tr><th scope=row>234261_at</th><td>-0.7375793</td><td> 4.591648</td><td>-5.617463</td><td>1.976787e-04</td><td>0.8484434</td><td>-3.087075</td></tr>
	<tr><th scope=row>232363_at</th><td> 0.7204464</td><td> 5.589081</td><td> 5.609625</td><td>1.998920e-04</td><td>0.8484434</td><td>-3.088461</td></tr>
	<tr><th scope=row>201324_at</th><td> 3.1314145</td><td> 7.011367</td><td> 5.485937</td><td>2.385607e-04</td><td>0.8484434</td><td>-3.110812</td></tr>
	<tr><th scope=row>1562249_at</th><td> 0.7558017</td><td> 3.635280</td><td> 5.483832</td><td>2.392842e-04</td><td>0.8484434</td><td>-3.111200</td></tr>
	<tr><th scope=row>202445_s_at</th><td> 0.7651195</td><td> 6.774628</td><td> 5.410852</td><td>2.658764e-04</td><td>0.8484434</td><td>-3.124829</td></tr>
	<tr><th scope=row>218559_s_at</th><td> 2.6383287</td><td>10.224166</td><td> 5.293459</td><td>3.154910e-04</td><td>0.8484434</td><td>-3.147449</td></tr>
	<tr><th scope=row>203543_s_at</th><td> 2.3900498</td><td> 9.127486</td><td> 5.256977</td><td>3.328537e-04</td><td>0.8484434</td><td>-3.154658</td></tr>
	<tr><th scope=row>222670_s_at</th><td> 1.8936729</td><td> 9.504714</td><td> 5.173428</td><td>3.765705e-04</td><td>0.8484434</td><td>-3.171497</td></tr>
	<tr><th scope=row>226974_at</th><td> 1.1430556</td><td> 7.216918</td><td> 5.157379</td><td>3.856479e-04</td><td>0.8484434</td><td>-3.174785</td></tr>
	<tr><th scope=row>233383_at</th><td>-0.4495825</td><td> 4.858983</td><td>-5.155069</td><td>3.869733e-04</td><td>0.8484434</td><td>-3.175259</td></tr>
	<tr><th scope=row>244353_s_at</th><td> 1.3289572</td><td> 6.833695</td><td> 5.150799</td><td>3.894365e-04</td><td>0.8484434</td><td>-3.176138</td></tr>
	<tr><th scope=row>1553391_at</th><td>-0.5431415</td><td> 4.139402</td><td>-5.099880</td><td>4.201297e-04</td><td>0.8484434</td><td>-3.186708</td></tr>
	<tr><th scope=row>231270_at</th><td> 0.4651273</td><td> 5.693016</td><td> 5.056856</td><td>4.480713e-04</td><td>0.8484434</td><td>-3.195778</td></tr>
	<tr><th scope=row>205841_at</th><td> 1.1069886</td><td> 7.759010</td><td> 5.035117</td><td>4.629356e-04</td><td>0.8484434</td><td>-3.200409</td></tr>
	<tr><th scope=row>230151_at</th><td> 0.6930578</td><td> 7.104013</td><td> 5.009513</td><td>4.811181e-04</td><td>0.8484434</td><td>-3.205907</td></tr>
	<tr><th scope=row>203962_s_at</th><td> 0.6734192</td><td>10.593867</td><td> 4.901058</td><td>5.669861e-04</td><td>0.8484434</td><td>-3.229713</td></tr>
	<tr><th scope=row>223169_s_at</th><td> 1.1723316</td><td> 7.710865</td><td> 4.872767</td><td>5.919663e-04</td><td>0.8484434</td><td>-3.236064</td></tr>
	<tr><th scope=row>227410_at</th><td> 1.8137747</td><td> 8.119707</td><td> 4.862006</td><td>6.017725e-04</td><td>0.8484434</td><td>-3.238495</td></tr>
	<tr><th scope=row>227569_at</th><td> 0.5997282</td><td>10.893658</td><td> 4.850481</td><td>6.124668e-04</td><td>0.8484434</td><td>-3.241108</td></tr>
	<tr><th scope=row>222877_at</th><td> 1.1503280</td><td> 4.873054</td><td> 4.798491</td><td>6.632744e-04</td><td>0.8484434</td><td>-3.253018</td></tr>
	<tr><th scope=row>1553284_s_at</th><td>-0.4645544</td><td> 5.170025</td><td>-4.778128</td><td>6.843758e-04</td><td>0.8484434</td><td>-3.257739</td></tr>
	<tr><th scope=row>227949_at</th><td> 0.7689198</td><td> 5.133790</td><td> 4.773179</td><td>6.896111e-04</td><td>0.8484434</td><td>-3.258891</td></tr>
	<tr><th scope=row>...</th><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
	<tr><th scope=row>215283_at</th><td>-2.403409e-04</td><td> 5.753796</td><td>-1.047986e-03</td><td>0.9991838</td><td>0.9997141</td><td>-4.809012</td></tr>
	<tr><th scope=row>227460_at</th><td>-1.534701e-04</td><td> 3.640915</td><td>-9.966703e-04</td><td>0.9992238</td><td>0.9997358</td><td>-4.809012</td></tr>
	<tr><th scope=row>215170_s_at</th><td>-3.292286e-04</td><td> 5.363076</td><td>-9.308910e-04</td><td>0.9992750</td><td>0.9997559</td><td>-4.809012</td></tr>
	<tr><th scope=row>208426_x_at</th><td> 7.893967e-05</td><td> 3.929755</td><td> 9.238792e-04</td><td>0.9992805</td><td>0.9997559</td><td>-4.809012</td></tr>
	<tr><th scope=row>201101_s_at</th><td> 1.418590e-04</td><td>11.118293</td><td> 8.786085e-04</td><td>0.9993157</td><td>0.9997729</td><td>-4.809012</td></tr>
	<tr><th scope=row>216005_at</th><td>-6.758614e-05</td><td> 3.894155</td><td>-8.485050e-04</td><td>0.9993392</td><td>0.9997778</td><td>-4.809012</td></tr>
	<tr><th scope=row>240890_at</th><td> 6.336304e-05</td><td> 3.871573</td><td> 8.191105e-04</td><td>0.9993621</td><td>0.9997778</td><td>-4.809012</td></tr>
	<tr><th scope=row>1556987_s_at</th><td>-1.195134e-04</td><td> 5.591520</td><td>-7.990722e-04</td><td>0.9993777</td><td>0.9997778</td><td>-4.809012</td></tr>
	<tr><th scope=row>1560525_at</th><td> 9.080890e-05</td><td> 4.479256</td><td> 7.783917e-04</td><td>0.9993938</td><td>0.9997778</td><td>-4.809012</td></tr>
	<tr><th scope=row>215117_at</th><td>-6.032394e-05</td><td> 3.210949</td><td>-6.521922e-04</td><td>0.9994921</td><td>0.9998430</td><td>-4.809012</td></tr>
	<tr><th scope=row>1555448_at</th><td>-1.302898e-04</td><td> 6.417103</td><td>-6.477286e-04</td><td>0.9994955</td><td>0.9998430</td><td>-4.809012</td></tr>
	<tr><th scope=row>1552586_at</th><td> 1.091079e-04</td><td> 5.570028</td><td> 5.512767e-04</td><td>0.9995707</td><td>0.9998527</td><td>-4.809013</td></tr>
	<tr><th scope=row>1566033_at</th><td>-3.930598e-05</td><td> 3.910951</td><td>-5.288197e-04</td><td>0.9995881</td><td>0.9998527</td><td>-4.809013</td></tr>
	<tr><th scope=row>1557098_s_at</th><td>-1.257625e-04</td><td> 6.062424</td><td>-5.174717e-04</td><td>0.9995970</td><td>0.9998527</td><td>-4.809013</td></tr>
	<tr><th scope=row>209812_x_at</th><td> 8.085186e-05</td><td> 6.773522</td><td> 5.069707e-04</td><td>0.9996052</td><td>0.9998527</td><td>-4.809013</td></tr>
	<tr><th scope=row>217027_x_at</th><td>-4.676716e-05</td><td> 9.199654</td><td>-4.817010e-04</td><td>0.9996248</td><td>0.9998527</td><td>-4.809013</td></tr>
	<tr><th scope=row>223427_s_at</th><td>-8.826246e-05</td><td> 7.532489</td><td>-4.668225e-04</td><td>0.9996364</td><td>0.9998527</td><td>-4.809013</td></tr>
	<tr><th scope=row>215666_at</th><td>-4.745100e-05</td><td> 3.513468</td><td>-4.494479e-04</td><td>0.9996500</td><td>0.9998527</td><td>-4.809013</td></tr>
	<tr><th scope=row>1568663_a_at</th><td> 5.849318e-05</td><td> 3.784688</td><td> 4.271066e-04</td><td>0.9996674</td><td>0.9998527</td><td>-4.809013</td></tr>
	<tr><th scope=row>235179_at</th><td> 9.864855e-05</td><td> 6.183387</td><td> 4.239390e-04</td><td>0.9996698</td><td>0.9998527</td><td>-4.809013</td></tr>
	<tr><th scope=row>208074_s_at</th><td> 5.595745e-05</td><td>11.862936</td><td> 3.286784e-04</td><td>0.9997440</td><td>0.9999086</td><td>-4.809013</td></tr>
	<tr><th scope=row>1554803_s_at</th><td> 4.574770e-05</td><td> 5.060954</td><td> 2.455111e-04</td><td>0.9998088</td><td>0.9999490</td><td>-4.809013</td></tr>
	<tr><th scope=row>1559228_at</th><td> 3.194653e-05</td><td> 4.166964</td><td> 2.230125e-04</td><td>0.9998263</td><td>0.9999490</td><td>-4.809013</td></tr>
	<tr><th scope=row>229217_at</th><td>-2.384651e-05</td><td> 6.346017</td><td>-1.740644e-04</td><td>0.9998644</td><td>0.9999490</td><td>-4.809013</td></tr>
	<tr><th scope=row>233856_at</th><td> 1.446460e-05</td><td> 3.595306</td><td> 1.628904e-04</td><td>0.9998731</td><td>0.9999490</td><td>-4.809013</td></tr>
	<tr><th scope=row>1566266_at</th><td>-1.372198e-05</td><td> 3.376231</td><td>-1.451925e-04</td><td>0.9998869</td><td>0.9999490</td><td>-4.809013</td></tr>
	<tr><th scope=row>216047_x_at</th><td>-1.793429e-05</td><td> 5.717270</td><td>-1.358965e-04</td><td>0.9998942</td><td>0.9999490</td><td>-4.809013</td></tr>
	<tr><th scope=row>226061_s_at</th><td>-6.476086e-06</td><td> 5.179055</td><td>-6.441839e-05</td><td>0.9999498</td><td>0.9999864</td><td>-4.809013</td></tr>
	<tr><th scope=row>229863_s_at</th><td> 1.966119e-06</td><td> 8.510457</td><td> 2.052237e-05</td><td>0.9999840</td><td>0.9999894</td><td>-4.809013</td></tr>
	<tr><th scope=row>243796_at</th><td>-3.021569e-06</td><td> 7.134296</td><td>-1.366763e-05</td><td>0.9999894</td><td>0.9999894</td><td>-4.809013</td></tr>
</tbody>
</table>




```R

```


```R

```


```R

```


```R

```


```R

```


```R

```
