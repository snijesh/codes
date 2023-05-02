
If you want to each file one by one use below codes

```
rsync -avP path/to/harddrive/ path/to/destination/
```

for example:

```
rsync -avP /media/home/hdd/R3278_R1.fastq.gz /media/home/edrive/data/R3278_R1.fastq.gz
rsync -avP /media/home/hdd/R3278_R2.fastq.gz /media/home/edrive/data/R3278_R2.fastq.gz
rsync -avP /media/home/hdd/R3279_R1.fastq.gz /media/home/edrive/data/R3279_R1.fastq.gz
```

If you want to copy all the files in the directory to another directory follow:
```
if rsync -avP /media/molmed/Data/sfari/samples/S1/ /media/molmed/sfaridata/fastq; then echo "transfer complete" else echo "transfer failed" fi
```




In shell script copy.sh
```
#!/bin/bash
if rsync -avP /media/molmed/Data/sfari/samples/S1/ /media/molmed/sfaridata/fastq; then
  echo "transfer complete"
else
  echo "transfer failed"
fi

```
