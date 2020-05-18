#!/bin/tcsh

set MAX_THREADS = 16

foreach threads (`seq ${MAX_THREADS}`)
    setenv OMP_NUM_THREADS $threads
    ./error
end
