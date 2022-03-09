# CGM-based Parallel Solutions For The MPP

This repository stores coarse-grained multicomputer parallel algorithms to solve the minimum cost parenthesization problem (MPP) using : 

1.  the regular partitioning techniques : 
    - Kechid and Myoupo [hal-01007646](https://hal.archives-ouvertes.fr/hal-01007646)
    - Kengne and Myoupo [doi:10.1007/s11227-011-0601-9](https://doi.org/10.1007/s11227-011-0601-9) 

2.  the irregular partitioning technique : Lacmou and Kengne [hal-01900171](https://hal.archives-ouvertes.fr/hal-01900171)

3.  the *k*-block splitting technique : Lacmou et *al.* [doi:10.1007/s11227-021-04069-9](https://doi.org/10.1007/s11227-021-04069-9) 

4.  the four-splitting technique.

## Build projects

Just execute in each root folder:

````
make
````

## Run projects

#### A. Configuration

The config file is located at the path :

````
{project}/resources/config/config.json
````

You can set the data set, the data size, the algorithm be used, etc.

````
...
...
"datasets": {
    "datasets-dir": "resources/datasets/",
    "datasets-type": "dimensions/",
    "datasets-dirname": "65537/",
    "datasets-filename": "65537",
    "datasets-extension": ".data",
    "datasets-max-number-matrix": "65537",
    "datasets-input": {
        "datasets-input-id": "_2",
        "datasets-input-max-fragmentation": "3",
        "datasets-input-max-number-matrix": "1024",
        "datasets-input-algorithm": "5",
        "datasets-input-partitioning-technique": "1"
    },
    "datasets-create": {
        "datasets-create-id": "_4"
    }
},
...
...
````

The algorithms be used are identifyied by a number :

0.  Godbole
1.  *k*-block splitting technique (*diagonal by diagonal* evaluation approach)
2.  *k*-block splitting technique (*k-block by k-block* evaluation approach)
3.  regular partitioning technique (Kechid and Myoupo)
4.  irregular partitioning technique (Lacmou and Kengne)
5.  four-splitting technique (Lacmou et *al.*)

The partitioning techniques (for the parallel solutions) are also identifyied by a number :

0.  the DAG is divided by *p* (Kechid and Myoupo)
1.  the DAG is divided by *sqrt(p)* (Kengne and Myoupo, Lacmou and Kengne, and Lacmou et *al.*)

#### B. Parallel algorithms

Run :

````
mpirun -np {number of processors} -hostfile hosts ./bin/CGM-MPP.run {data size} {number of fragmentations} {algorithm be used} {partitioning technique}
````

For instance :

````
mpirun -np 8 -hostfile hosts ./bin/CGM-MPP.run 10240 1 5 1
````

## License
[MIT](LICENSE)