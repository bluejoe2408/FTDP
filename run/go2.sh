#!/bin/sh
#PJM -N "test2"
#PJM -L rscgrp=lecture
#PJM -L node=1
#PJM --omp thread=48
#PJM -L elapse=00:15:00
#PJM -g gt69
#PJM -j
#PJM -e err
#PJM -o fin2f1.lst

./sol20
./sol20
./sol20
./sol20
./sol20
