#!/bin/bash

for pf1 in 0.01 0.05 0.01 0.2
do
        for er in 5 10 15 20 25 35 45 50
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
