#!/bin/sh
#PJM -N "test2"
#PJM -L rscgrp=lecture
#PJM -L node=1
#PJM --omp thread=24
#PJM -L elapse=00:15:00
#PJM -g gt00
#PJM -j
#PJM -e err
#PJM -o test4.lst

export KMP_AFFINITY=verbose
./sol20

export KMP_AFFINITY=granularity=fine,compact,verbose
./sol20

export KMP_AFFINITY=granularity=fine,balanced,verbose
./sol20

export KMP_AFFINITY=granularity=fine,scatter,verbose
./sol20
