#!/bin/bash

for pf1 in 0.05 0.1 0.2
do
        for er in 50 40 20 30 10 5
        do
                for delta in 0.1 0.2
                do
                        for beta in 2 4
                        do
                                sbatch submit_job.sh $beta $pf1 $delta $er
                        done
                done
        done
done
