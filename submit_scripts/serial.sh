#!/bin/bash

# declare a name for this job (replace <jobname> with something more descriptive)
#$ -N Schultz_Lane

# request the queue for this job
# replace <queue_name> with queue_name e.g class.q
#$ -q class.q

# request 4 hours of wall time
#$ -l h_rt=04:00:00

# run the job from the directory of submission. Uncomment only if you don't want the defaults.
#$ -cwd
# combine SGE standard output and error files
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
# transfer all your environment variables. Uncomment only if you don't want the defults
#$ -V

./a.out

