# BEACON v0.1
BEACON is a novel method to query the relatedness of a candidate gene network by measuring the degrees of separation of protein-protein functional associations from the STRING database


## Run Guide 
Run guide:
1. First, please download all protein links from STRING. Ex:
```wget \\"https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz\\" ```                                                          
2. Next, run beaconPathfind.py to generate a cache of the shortest number of functional links between each pair of proteins in the database. Ex: 
```python beaconPathfind.py --input 9606.protein.links.v12.0.txt.gz --threshold 400 --output cache```                                                                
3. Generate a BEACON score and compare it to random protein lists using run_BEACON.pl. Ex:                                                                     
```perl run_BEACON.pl --input proteins_of_interest.txt --cache cache_directory --output output_directory\'94```


## Usage

```USAGE: perl run_BEACON.pl proteins_of_interest.txt cache_directory  ```}
