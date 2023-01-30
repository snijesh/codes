```python
import os
import pandas as pd
import glob
```


```python
os.chdir('ashtest/')
```


```python
files = glob.glob('*.txt')
files
```




    ['107.txt', '78.txt']




```python
for i in files:
    data=pd.read_table(i, sep="\t", low_memory=False)
    data=data.assign(SampleID=i.replace('.txt', ''))
    data.to_csv(i, sep="\t", index=0)
```

    C:\Users\snijesh\AppData\Local\Temp\ipykernel_10296\4194162390.py:2: DtypeWarning: Columns (253,254,258) have mixed types. Specify dtype option on import or set low_memory=False.
      data=pd.read_table(i, sep="\t")
    C:\Users\snijesh\AppData\Local\Temp\ipykernel_10296\4194162390.py:2: DtypeWarning: Columns (253,254,258) have mixed types. Specify dtype option on import or set low_memory=False.
      data=pd.read_table(i, sep="\t")
    


```python
!type *.txt > combined.txt ## use this in windows
```

    
    107.txt
    
    
    
    78.txt
    
    
    
    combined.txt
    
    
    


```python
# !cat *.txt > combined.txt 
# use this in linux
```


```python
#read the combined vcf file
vcf = pd.read_table('combined.txt', sep='\t', low_memory=False)
vcf
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Chr</th>
      <th>Start</th>
      <th>End</th>
      <th>Ref</th>
      <th>Alt</th>
      <th>Func.refGene</th>
      <th>Gene.refGene</th>
      <th>GeneDetail.refGene</th>
      <th>ExonicFunc.refGene</th>
      <th>AAChange.refGene</th>
      <th>...</th>
      <th>Otherinfo5</th>
      <th>Otherinfo6</th>
      <th>Otherinfo7</th>
      <th>Otherinfo8</th>
      <th>Otherinfo9</th>
      <th>Otherinfo10</th>
      <th>Otherinfo11</th>
      <th>Otherinfo12</th>
      <th>Otherinfo13</th>
      <th>SampleID</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>chr1</td>
      <td>12198</td>
      <td>12198</td>
      <td>G</td>
      <td>C</td>
      <td>ncRNA_exonic</td>
      <td>DDX11L1;LOC102725121</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>12198</td>
      <td>.</td>
      <td>G</td>
      <td>C</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=36;WT=0;HET=0;HOM=1;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>1/1:139:37:36:6:30:83.33%:1.1855E-14:35:41:2:4...</td>
      <td>107</td>
    </tr>
    <tr>
      <th>1</th>
      <td>chr1</td>
      <td>12783</td>
      <td>12783</td>
      <td>G</td>
      <td>A</td>
      <td>ncRNA_intronic</td>
      <td>DDX11L1;LOC102725121</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>12783</td>
      <td>.</td>
      <td>G</td>
      <td>A</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=22;WT=0;HET=1;HOM=0;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>0/1:45:22:22:10:12:54.55%:3.066E-5:41:30:6:4:8:4</td>
      <td>107</td>
    </tr>
    <tr>
      <th>2</th>
      <td>chr1</td>
      <td>13302</td>
      <td>13302</td>
      <td>C</td>
      <td>T</td>
      <td>ncRNA_exonic</td>
      <td>DDX11L1;LOC102725121</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>13302</td>
      <td>.</td>
      <td>C</td>
      <td>T</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=15;WT=0;HET=1;HOM=0;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>0/1:46:15:15:4:11:73.33%:2.4988E-5:33:40:3:1:9:2</td>
      <td>107</td>
    </tr>
    <tr>
      <th>3</th>
      <td>chr1</td>
      <td>13896</td>
      <td>13896</td>
      <td>C</td>
      <td>A</td>
      <td>ncRNA_exonic</td>
      <td>DDX11L1;LOC102725121</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>13896</td>
      <td>.</td>
      <td>C</td>
      <td>A</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=31;WT=0;HET=1;HOM=0;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>0/1:37:31:31:20:11:35.48%:1.6659E-4:38:38:15:5...</td>
      <td>107</td>
    </tr>
    <tr>
      <th>4</th>
      <td>chr1</td>
      <td>14522</td>
      <td>14522</td>
      <td>G</td>
      <td>A</td>
      <td>ncRNA_exonic</td>
      <td>WASH7P</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>14522</td>
      <td>.</td>
      <td>G</td>
      <td>A</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=66;WT=0;HET=1;HOM=0;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>0/1:67:70:66:46:20:30.3%:1.7467E-7:38:33:22:24...</td>
      <td>107</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>230617</th>
      <td>chrM</td>
      <td>15933</td>
      <td>15933</td>
      <td>C</td>
      <td>T</td>
      <td>intergenic</td>
      <td>MIR12136;NONE</td>
      <td>dist=8418;dist=NONE</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>15933</td>
      <td>.</td>
      <td>C</td>
      <td>T</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=21;WT=0;HET=0;HOM=1;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>1/1:117:21:21:0:21:100%:1.8578E-12:0:39:0:0:17:4</td>
      <td>78</td>
    </tr>
    <tr>
      <th>230618</th>
      <td>chrM</td>
      <td>16173</td>
      <td>16173</td>
      <td>C</td>
      <td>T</td>
      <td>intergenic</td>
      <td>MIR12136;NONE</td>
      <td>dist=8658;dist=NONE</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>16173</td>
      <td>.</td>
      <td>C</td>
      <td>T</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=11;WT=0;HET=0;HOM=1;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>1/1:58:11:11:0:11:100%:1.4176E-6:0:42:0:0:7:4</td>
      <td>78</td>
    </tr>
    <tr>
      <th>230619</th>
      <td>chrM</td>
      <td>16223</td>
      <td>16223</td>
      <td>C</td>
      <td>T</td>
      <td>intergenic</td>
      <td>MIR12136;NONE</td>
      <td>dist=8708;dist=NONE</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>16223</td>
      <td>.</td>
      <td>C</td>
      <td>T</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=10;WT=0;HET=0;HOM=1;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>1/1:52:10:10:0:10:100%:5.4125E-6:0:50:0:0:4:6</td>
      <td>78</td>
    </tr>
    <tr>
      <th>230620</th>
      <td>chrM</td>
      <td>16225</td>
      <td>16225</td>
      <td>T</td>
      <td>C</td>
      <td>intergenic</td>
      <td>MIR12136;NONE</td>
      <td>dist=8710;dist=NONE</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>16225</td>
      <td>.</td>
      <td>T</td>
      <td>C</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=10;WT=0;HET=0;HOM=1;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>1/1:52:10:10:0:10:100%:5.4125E-6:0:52:0:0:4:6</td>
      <td>78</td>
    </tr>
    <tr>
      <th>230621</th>
      <td>Chr</td>
      <td>Start</td>
      <td>End</td>
      <td>Ref</td>
      <td>Alt</td>
      <td>Func.refGene</td>
      <td>Gene.refGen</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
<p>230622 rows × 273 columns</p>
</div>




```python
#remove the column names of multiple files
vcf=vcf[vcf['Chr']!='Chr']
vcf
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Chr</th>
      <th>Start</th>
      <th>End</th>
      <th>Ref</th>
      <th>Alt</th>
      <th>Func.refGene</th>
      <th>Gene.refGene</th>
      <th>GeneDetail.refGene</th>
      <th>ExonicFunc.refGene</th>
      <th>AAChange.refGene</th>
      <th>...</th>
      <th>Otherinfo5</th>
      <th>Otherinfo6</th>
      <th>Otherinfo7</th>
      <th>Otherinfo8</th>
      <th>Otherinfo9</th>
      <th>Otherinfo10</th>
      <th>Otherinfo11</th>
      <th>Otherinfo12</th>
      <th>Otherinfo13</th>
      <th>SampleID</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>chr1</td>
      <td>12198</td>
      <td>12198</td>
      <td>G</td>
      <td>C</td>
      <td>ncRNA_exonic</td>
      <td>DDX11L1;LOC102725121</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>12198</td>
      <td>.</td>
      <td>G</td>
      <td>C</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=36;WT=0;HET=0;HOM=1;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>1/1:139:37:36:6:30:83.33%:1.1855E-14:35:41:2:4...</td>
      <td>107</td>
    </tr>
    <tr>
      <th>1</th>
      <td>chr1</td>
      <td>12783</td>
      <td>12783</td>
      <td>G</td>
      <td>A</td>
      <td>ncRNA_intronic</td>
      <td>DDX11L1;LOC102725121</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>12783</td>
      <td>.</td>
      <td>G</td>
      <td>A</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=22;WT=0;HET=1;HOM=0;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>0/1:45:22:22:10:12:54.55%:3.066E-5:41:30:6:4:8:4</td>
      <td>107</td>
    </tr>
    <tr>
      <th>2</th>
      <td>chr1</td>
      <td>13302</td>
      <td>13302</td>
      <td>C</td>
      <td>T</td>
      <td>ncRNA_exonic</td>
      <td>DDX11L1;LOC102725121</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>13302</td>
      <td>.</td>
      <td>C</td>
      <td>T</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=15;WT=0;HET=1;HOM=0;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>0/1:46:15:15:4:11:73.33%:2.4988E-5:33:40:3:1:9:2</td>
      <td>107</td>
    </tr>
    <tr>
      <th>3</th>
      <td>chr1</td>
      <td>13896</td>
      <td>13896</td>
      <td>C</td>
      <td>A</td>
      <td>ncRNA_exonic</td>
      <td>DDX11L1;LOC102725121</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>13896</td>
      <td>.</td>
      <td>C</td>
      <td>A</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=31;WT=0;HET=1;HOM=0;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>0/1:37:31:31:20:11:35.48%:1.6659E-4:38:38:15:5...</td>
      <td>107</td>
    </tr>
    <tr>
      <th>4</th>
      <td>chr1</td>
      <td>14522</td>
      <td>14522</td>
      <td>G</td>
      <td>A</td>
      <td>ncRNA_exonic</td>
      <td>WASH7P</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>14522</td>
      <td>.</td>
      <td>G</td>
      <td>A</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=66;WT=0;HET=1;HOM=0;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>0/1:67:70:66:46:20:30.3%:1.7467E-7:38:33:22:24...</td>
      <td>107</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>230616</th>
      <td>chrM</td>
      <td>15302</td>
      <td>15302</td>
      <td>A</td>
      <td>G</td>
      <td>intergenic</td>
      <td>MIR12136;NONE</td>
      <td>dist=7787;dist=NONE</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>15302</td>
      <td>.</td>
      <td>A</td>
      <td>G</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=14;WT=0;HET=0;HOM=1;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>1/1:76:14:14:0:14:100%:2.4927E-8:0:42:0:0:10:4</td>
      <td>78</td>
    </tr>
    <tr>
      <th>230617</th>
      <td>chrM</td>
      <td>15933</td>
      <td>15933</td>
      <td>C</td>
      <td>T</td>
      <td>intergenic</td>
      <td>MIR12136;NONE</td>
      <td>dist=8418;dist=NONE</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>15933</td>
      <td>.</td>
      <td>C</td>
      <td>T</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=21;WT=0;HET=0;HOM=1;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>1/1:117:21:21:0:21:100%:1.8578E-12:0:39:0:0:17:4</td>
      <td>78</td>
    </tr>
    <tr>
      <th>230618</th>
      <td>chrM</td>
      <td>16173</td>
      <td>16173</td>
      <td>C</td>
      <td>T</td>
      <td>intergenic</td>
      <td>MIR12136;NONE</td>
      <td>dist=8658;dist=NONE</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>16173</td>
      <td>.</td>
      <td>C</td>
      <td>T</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=11;WT=0;HET=0;HOM=1;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>1/1:58:11:11:0:11:100%:1.4176E-6:0:42:0:0:7:4</td>
      <td>78</td>
    </tr>
    <tr>
      <th>230619</th>
      <td>chrM</td>
      <td>16223</td>
      <td>16223</td>
      <td>C</td>
      <td>T</td>
      <td>intergenic</td>
      <td>MIR12136;NONE</td>
      <td>dist=8708;dist=NONE</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>16223</td>
      <td>.</td>
      <td>C</td>
      <td>T</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=10;WT=0;HET=0;HOM=1;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>1/1:52:10:10:0:10:100%:5.4125E-6:0:50:0:0:4:6</td>
      <td>78</td>
    </tr>
    <tr>
      <th>230620</th>
      <td>chrM</td>
      <td>16225</td>
      <td>16225</td>
      <td>T</td>
      <td>C</td>
      <td>intergenic</td>
      <td>MIR12136;NONE</td>
      <td>dist=8710;dist=NONE</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>16225</td>
      <td>.</td>
      <td>T</td>
      <td>C</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=10;WT=0;HET=0;HOM=1;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>1/1:52:10:10:0:10:100%:5.4125E-6:0:52:0:0:4:6</td>
      <td>78</td>
    </tr>
  </tbody>
</table>
<p>230618 rows × 273 columns</p>
</div>




```python

```


```python
# filter by gene of interest (contains match)
DDX11L1 = vcf[vcf['Gene.refGene'].str.contains("DDX11L1")]
DDX11L1
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Chr</th>
      <th>Start</th>
      <th>End</th>
      <th>Ref</th>
      <th>Alt</th>
      <th>Func.refGene</th>
      <th>Gene.refGene</th>
      <th>GeneDetail.refGene</th>
      <th>ExonicFunc.refGene</th>
      <th>AAChange.refGene</th>
      <th>...</th>
      <th>Otherinfo5</th>
      <th>Otherinfo6</th>
      <th>Otherinfo7</th>
      <th>Otherinfo8</th>
      <th>Otherinfo9</th>
      <th>Otherinfo10</th>
      <th>Otherinfo11</th>
      <th>Otherinfo12</th>
      <th>Otherinfo13</th>
      <th>SampleID</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>chr1</td>
      <td>12198</td>
      <td>12198</td>
      <td>G</td>
      <td>C</td>
      <td>ncRNA_exonic</td>
      <td>DDX11L1;LOC102725121</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>12198</td>
      <td>.</td>
      <td>G</td>
      <td>C</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=36;WT=0;HET=0;HOM=1;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>1/1:139:37:36:6:30:83.33%:1.1855E-14:35:41:2:4...</td>
      <td>107</td>
    </tr>
    <tr>
      <th>1</th>
      <td>chr1</td>
      <td>12783</td>
      <td>12783</td>
      <td>G</td>
      <td>A</td>
      <td>ncRNA_intronic</td>
      <td>DDX11L1;LOC102725121</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>12783</td>
      <td>.</td>
      <td>G</td>
      <td>A</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=22;WT=0;HET=1;HOM=0;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>0/1:45:22:22:10:12:54.55%:3.066E-5:41:30:6:4:8:4</td>
      <td>107</td>
    </tr>
    <tr>
      <th>2</th>
      <td>chr1</td>
      <td>13302</td>
      <td>13302</td>
      <td>C</td>
      <td>T</td>
      <td>ncRNA_exonic</td>
      <td>DDX11L1;LOC102725121</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>13302</td>
      <td>.</td>
      <td>C</td>
      <td>T</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=15;WT=0;HET=1;HOM=0;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>0/1:46:15:15:4:11:73.33%:2.4988E-5:33:40:3:1:9:2</td>
      <td>107</td>
    </tr>
    <tr>
      <th>3</th>
      <td>chr1</td>
      <td>13896</td>
      <td>13896</td>
      <td>C</td>
      <td>A</td>
      <td>ncRNA_exonic</td>
      <td>DDX11L1;LOC102725121</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>13896</td>
      <td>.</td>
      <td>C</td>
      <td>A</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=31;WT=0;HET=1;HOM=0;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>0/1:37:31:31:20:11:35.48%:1.6659E-4:38:38:15:5...</td>
      <td>107</td>
    </tr>
    <tr>
      <th>33789</th>
      <td>chr16</td>
      <td>61651</td>
      <td>61651</td>
      <td>A</td>
      <td>G</td>
      <td>ncRNA_exonic</td>
      <td>DDX11L10</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>61651</td>
      <td>.</td>
      <td>A</td>
      <td>G</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=10;WT=0;HET=0;HOM=1;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>1/1:52:11:10:0:10:100%:5.4125E-6:0:37:0:0:10:0</td>
      <td>107</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>230589</th>
      <td>chrY</td>
      <td>59355764</td>
      <td>59355764</td>
      <td>A</td>
      <td>G</td>
      <td>intergenic</td>
      <td>WASIR1;DDX11L16</td>
      <td>dist=6263;dist=2565</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>59355764</td>
      <td>.</td>
      <td>A</td>
      <td>G</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=10;WT=0;HET=0;HOM=1;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>1/1:42:10:10:1:9:90%:5.9538E-5:40:44:1:0:6:3</td>
      <td>78</td>
    </tr>
    <tr>
      <th>230590</th>
      <td>chrY</td>
      <td>59357018</td>
      <td>59357018</td>
      <td>T</td>
      <td>C</td>
      <td>intergenic</td>
      <td>WASIR1;DDX11L16</td>
      <td>dist=7517;dist=1311</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>59357018</td>
      <td>.</td>
      <td>T</td>
      <td>C</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=21;WT=0;HET=1;HOM=0;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>0/1:36:21:21:11:10:47.62%:2.3971E-4:48:38:5:6:1:9</td>
      <td>78</td>
    </tr>
    <tr>
      <th>230591</th>
      <td>chrY</td>
      <td>59357695</td>
      <td>59357695</td>
      <td>C</td>
      <td>T</td>
      <td>downstream</td>
      <td>DDX11L16</td>
      <td>dist=634</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>59357695</td>
      <td>.</td>
      <td>C</td>
      <td>T</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=66;WT=0;HET=1;HOM=0;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>0/1:60:66:66:48:18:27.27%:9.974E-7:44:37:29:19...</td>
      <td>78</td>
    </tr>
    <tr>
      <th>230592</th>
      <td>chrY</td>
      <td>59357787</td>
      <td>59357787</td>
      <td>G</td>
      <td>A</td>
      <td>downstream</td>
      <td>DDX11L16</td>
      <td>dist=542</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>59357787</td>
      <td>.</td>
      <td>G</td>
      <td>A</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=38;WT=0;HET=1;HOM=0;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>0/1:40:38:38:26:12:31.58%:8.7276E-5:36:41:12:1...</td>
      <td>78</td>
    </tr>
    <tr>
      <th>230593</th>
      <td>chrY</td>
      <td>59358392</td>
      <td>59358392</td>
      <td>G</td>
      <td>T</td>
      <td>ncRNA_exonic</td>
      <td>DDX11L16</td>
      <td>.</td>
      <td>.</td>
      <td>.</td>
      <td>...</td>
      <td>59358392</td>
      <td>.</td>
      <td>G</td>
      <td>T</td>
      <td>.</td>
      <td>PASS</td>
      <td>ADP=38;WT=0;HET=1;HOM=0;NC=0</td>
      <td>GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:A...</td>
      <td>0/1:36:38:38:27:11:28.95%:2.1011E-4:44:39:12:1...</td>
      <td>78</td>
    </tr>
  </tbody>
</table>
<p>96 rows × 273 columns</p>
</div>




```python
#rfilter by gene of interest (exact match)
vcf=vcf[vcf['Gene.refGene']=="DDX11L1]
vcf
```


```python

```


```python

```


```python

```
