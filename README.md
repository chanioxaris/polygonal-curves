## Overview
This is an implementation of clustering algorithms for polygonal curves. User can choose between the below methods:

#### Initialization
1. K-means++
2. Random selection of k-points - (Simplest)
#### Assignment
1. Lloyd’s assignment - (Simplest)
2. Assignment by Range search (LSH)
#### Update
1. Mean Discrete Frechet curve
2. Partitioning Around Medoids (PAM) – (Improved)

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

## Output data files

### Output file 
The format of output text file, described by the following structure:
```
Algorithm: Ι<x>A<y>U<z>
Metric: Frechet or DTW
CLUSTER-1 {size: <int>, centroid: <curve_id>}
...
CLUSTER-k {size: <int>, centroid: <curve_id>}
clustering_time: <double> 
Silhouette: [s1,...,sk,stotal]
```
where ```x``` the chosen Initialization method (1 or 2), ```y``` the chosen Assignment method (1 or 2), ```z``` the chosen Update method (1 or 2)

### Output file (-complete)
If the user choose ```-complete``` flag then more information become available about the clusters:
```
CLUSTER-1 {curve_idA, curve_idB, ..., curve_idC}
CLUSTER-2 {curve_idA, curve_idB, ..., curve_idC}
...
CLUSTER-k {curve_idR, curve_idT, ..., curve_idZ}
```

### Evaluation (Silhouette)
Silhouette refers to a method of interpretation and validation of consistency within clusters of data. The technique provides a succinct graphical representation of how well each object lies within its cluster.
The silhouette value is a measure of how similar an object is to its own cluster (cohesion) compared to other clusters (separation). The silhouette ranges from −1 to +1, where a high value indicates that the object is well matched to its own cluster and poorly matched to neighboring clusters. If most objects have a high value, then the clustering configuration is appropriate. If many points have a low or negative value, then the clustering configuration may have too many or too few clusters.


![Silhouette](https://github.com/chanioxaris/Kmeans-Kmedois-PolygonalCurves/blob/master/img/silhouette.jpg)


## Compile

`./makefile`

## Usage

`./clustering –i [input file] –c [configuration file] -ο [output file] –d {Frechet, DTW} -complete (optional)`
