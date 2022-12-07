```python
import pandas as pd
import numpy as np
```


```python
#read the input data
df = pd.read_table('complex/input.txt', sep='\t')
#print number of columns and rows
print(df.shape)
#view first few lines of the dataframe
df.head()
```

    (1087, 17)
    




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
      <th>Sample</th>
      <th>Dnumber</th>
      <th>FamilyNumber</th>
      <th>DiagnosisRev</th>
      <th>chr_position</th>
      <th>Gene</th>
      <th>ExonicFunc</th>
      <th>avsnp150</th>
      <th>AF_sas</th>
      <th>AF_afr</th>
      <th>AF_amr</th>
      <th>AF_eas</th>
      <th>variant_class</th>
      <th>SFARI_gene_score</th>
      <th>syndromic</th>
      <th>genetic_category</th>
      <th>eagle</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>4538</td>
      <td>D638</td>
      <td>275</td>
      <td>Familial_control</td>
      <td>chr16_8866767_8866767_C_T</td>
      <td>ABAT</td>
      <td>nonsynonymous SNV</td>
      <td>rs555183419</td>
      <td>0.000200</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>missense</td>
      <td>2</td>
      <td>0</td>
      <td>Rare Single Gene Mutation, Genetic Association</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>14894</td>
      <td>D336</td>
      <td>9773</td>
      <td>Schizophrenia</td>
      <td>chr7_48567884_48567884_T_G</td>
      <td>ABCA13</td>
      <td>nonsynonymous SNV</td>
      <td>rs749109867</td>
      <td>0.000300</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>missense</td>
      <td>2</td>
      <td>0</td>
      <td>Rare Single Gene Mutation, Functional</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>16610</td>
      <td>D948</td>
      <td>9829</td>
      <td>Familial_control</td>
      <td>chr7_48412006_48412006_C_-</td>
      <td>ABCA13</td>
      <td>frameshift deletion</td>
      <td>NotAvailable</td>
      <td>0.000033</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>indel</td>
      <td>2</td>
      <td>0</td>
      <td>Rare Single Gene Mutation, Functional</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>14970</td>
      <td>D223</td>
      <td>9829</td>
      <td>Schizophrenia</td>
      <td>chr7_48412006_48412006_C_-</td>
      <td>ABCA13</td>
      <td>frameshift deletion</td>
      <td>NotAvailable</td>
      <td>0.000033</td>
      <td>0.000000</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>indel</td>
      <td>2</td>
      <td>0</td>
      <td>Rare Single Gene Mutation, Functional</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>15626</td>
      <td>D300</td>
      <td>10198</td>
      <td>Schizophrenia</td>
      <td>chr7_48411788_48411788_T_A</td>
      <td>ABCA13</td>
      <td>nonsynonymous SNV</td>
      <td>rs771610168</td>
      <td>0.000000</td>
      <td>0.000065</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>missense</td>
      <td>2</td>
      <td>0</td>
      <td>Rare Single Gene Mutation, Functional</td>
      <td>0.0</td>
    </tr>
  </tbody>
</table>
</div>




```python
#to see the column names
df.columns
```




    Index(['Sample', 'Dnumber', 'FamilyNumber', 'DiagnosisRev', 'chr_position',
           'Gene', 'ExonicFunc', 'avsnp150', 'AF_sas', 'AF_afr', 'AF_amr',
           'AF_eas', 'variant_class', 'SFARI_gene_score', 'syndromic',
           'genetic_category', 'eagle'],
          dtype='object')




```python
#shrink data of interest
df1 = df.groupby(['Gene','chr_position', 'FamilyNumber'])[['Sample','Dnumber','ExonicFunc','avsnp150','variant_class', 'genetic_category']].agg(lambda x: '|'.join(x.unique()))
df1.head()
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
      <th></th>
      <th></th>
      <th>Sample</th>
      <th>Dnumber</th>
      <th>ExonicFunc</th>
      <th>avsnp150</th>
      <th>variant_class</th>
      <th>genetic_category</th>
    </tr>
    <tr>
      <th>Gene</th>
      <th>chr_position</th>
      <th>FamilyNumber</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ABAT</th>
      <th>chr16_8866767_8866767_C_T</th>
      <th>275</th>
      <td>4538</td>
      <td>D638</td>
      <td>nonsynonymous SNV</td>
      <td>rs555183419</td>
      <td>missense</td>
      <td>Rare Single Gene Mutation, Genetic Association</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">ABCA13</th>
      <th>chr7_48411788_48411788_T_A</th>
      <th>10198</th>
      <td>15626</td>
      <td>D300</td>
      <td>nonsynonymous SNV</td>
      <td>rs771610168</td>
      <td>missense</td>
      <td>Rare Single Gene Mutation, Functional</td>
    </tr>
    <tr>
      <th>chr7_48412006_48412006_C_-</th>
      <th>9829</th>
      <td>16610|14970</td>
      <td>D948|D223</td>
      <td>frameshift deletion</td>
      <td>NotAvailable</td>
      <td>indel</td>
      <td>Rare Single Gene Mutation, Functional</td>
    </tr>
    <tr>
      <th>chr7_48567884_48567884_T_G</th>
      <th>9773</th>
      <td>14894</td>
      <td>D336</td>
      <td>nonsynonymous SNV</td>
      <td>rs749109867</td>
      <td>missense</td>
      <td>Rare Single Gene Mutation, Functional</td>
    </tr>
    <tr>
      <th>ABCA7</th>
      <th>chr19_1055224_1055224_C_T</th>
      <th>gg30</th>
      <td>ICNBSX329</td>
      <td>NotAvailable</td>
      <td>nonsynonymous SNV</td>
      <td>NotAvailable</td>
      <td>missense</td>
      <td>Rare Single Gene Mutation</td>
    </tr>
  </tbody>
</table>
</div>




```python
#create count for DiagnosisRev
df2 = pd.get_dummies(df.set_index(['Gene','chr_position', 'FamilyNumber'])['DiagnosisRev']).sum(level=[0, 1, 2])
df2.head()
```

    C:\Users\snijesh\AppData\Local\Temp\ipykernel_18256\120868633.py:2: FutureWarning: Using the level keyword in DataFrame and Series aggregations is deprecated and will be removed in a future version. Use groupby instead. df.sum(level=1) should use df.groupby(level=1).sum().
      df2 = pd.get_dummies(df.set_index(['Gene','chr_position', 'FamilyNumber'])['DiagnosisRev']).sum(level=[0, 1, 2])
    




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
      <th></th>
      <th></th>
      <th>Addiction</th>
      <th>BPAD</th>
      <th>Dementia</th>
      <th>Depression</th>
      <th>Familial_control</th>
      <th>OCD</th>
      <th>Schizophrenia</th>
      <th>population_control</th>
    </tr>
    <tr>
      <th>Gene</th>
      <th>chr_position</th>
      <th>FamilyNumber</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ABAT</th>
      <th>chr16_8866767_8866767_C_T</th>
      <th>275</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">ABCA13</th>
      <th>chr7_48567884_48567884_T_G</th>
      <th>9773</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>chr7_48412006_48412006_C_-</th>
      <th>9829</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>chr7_48411788_48411788_T_A</th>
      <th>10198</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ABCA7</th>
      <th>chr19_1058262_1058262_C_T</th>
      <th>7862</th>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
</div>




```python
#merge shrinked table and count table
df3 = pd.concat([df1, df2], axis=1)
df3.head()
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
      <th></th>
      <th></th>
      <th>Sample</th>
      <th>Dnumber</th>
      <th>ExonicFunc</th>
      <th>avsnp150</th>
      <th>variant_class</th>
      <th>genetic_category</th>
      <th>Addiction</th>
      <th>BPAD</th>
      <th>Dementia</th>
      <th>Depression</th>
      <th>Familial_control</th>
      <th>OCD</th>
      <th>Schizophrenia</th>
      <th>population_control</th>
    </tr>
    <tr>
      <th>Gene</th>
      <th>chr_position</th>
      <th>FamilyNumber</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>ABAT</th>
      <th>chr16_8866767_8866767_C_T</th>
      <th>275</th>
      <td>4538</td>
      <td>D638</td>
      <td>nonsynonymous SNV</td>
      <td>rs555183419</td>
      <td>missense</td>
      <td>Rare Single Gene Mutation, Genetic Association</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">ABCA13</th>
      <th>chr7_48411788_48411788_T_A</th>
      <th>10198</th>
      <td>15626</td>
      <td>D300</td>
      <td>nonsynonymous SNV</td>
      <td>rs771610168</td>
      <td>missense</td>
      <td>Rare Single Gene Mutation, Functional</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>chr7_48412006_48412006_C_-</th>
      <th>9829</th>
      <td>16610|14970</td>
      <td>D948|D223</td>
      <td>frameshift deletion</td>
      <td>NotAvailable</td>
      <td>indel</td>
      <td>Rare Single Gene Mutation, Functional</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>chr7_48567884_48567884_T_G</th>
      <th>9773</th>
      <td>14894</td>
      <td>D336</td>
      <td>nonsynonymous SNV</td>
      <td>rs749109867</td>
      <td>missense</td>
      <td>Rare Single Gene Mutation, Functional</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>ABCA7</th>
      <th>chr19_1055224_1055224_C_T</th>
      <th>gg30</th>
      <td>ICNBSX329</td>
      <td>NotAvailable</td>
      <td>nonsynonymous SNV</td>
      <td>NotAvailable</td>
      <td>missense</td>
      <td>Rare Single Gene Mutation</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>




```python
#write the output to a file
df3.to_csv('complex/output_final.txt', sep='\t')
```


```python

```


```python

```


```python

```


```python

```


```python

```


```python

```


```python

```


```python

```
