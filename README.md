## Overview


## Input data files

### Input file 
The format of input text file, described by the following structure:
```
@dimension {2,3,4} 
curve_id1 m1 (x11,y11) (x12,y12) (x13,y13) ... (x1m1,y1m1)
curve_id2 m2 (x21,y21) (x22,y22) (x23,y23) ... (x2m2,y2m2)
curve_id3 m3 (x31,y31) (x32,y32) (x33,y33) ... (x3m3,y3m3)
...
curve_idN mN (xN1,yN1) (xN2,yN2) (xN3,yN3) ... (xNmN,yNmN)
```
where ```mi``` the total points included in i curve, ```(xij,yij)``` the coordinates of point j in i curve, when dimesion is equal to 2

### Configuration file 
The format of configuration file, described by the following structure:
```
number_of_clusters: <int> 
number_of_grid_curves: <int> 
number_of_hash_tables: <int>
```

## Compile

`./makefile`

## Usage

`./clustering –i [input file] –c [configuration file] -ο [output file] –d {Frechet, DTW} -complete (optional)`
