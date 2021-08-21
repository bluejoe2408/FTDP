#!/bin/sh
#PJM -N "test2"
#PJM -L rscgrp=lecture
#PJM -L node=1
#PJM --omp thread=24
#PJM -L elapse=00:15:00
#PJM -g gt00
#PJM -j
#PJM -e err
#PJM -o test3.lst

./sol20
./sol20
./sol20
./sol20
./sol20

export KMP_AFFINITY=granularity=fine,compact
./sol20
./sol20
./sol20
./sol20
./sol20

export KMP_AFFINITY=granularity=fine,balanced
./sol20
./sol20
./sol20
./sol20
./sol20

export KMP_AFFINITY=granularity=fine,scatter
./sol20
./sol20
./sol20
./sol20
./sol20
