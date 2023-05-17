#!/usr/local/bin/bash -l

# after all jobs launched by vcfByGene.ipynb have finished, check if any encountered memory errors

JOB_IDS=$(dx find jobs --tag my_job_tag --brief -n 2000) # set n >= max jobs w/ tag

for job in $JOB_IDS; do
  # echo $job
  dx watch $job -l\
    ALERT --no-job-info --no-timestamps --color off -q |\
    grep -q "oom-killer"
  if [ "$?" == "0" ]; then
    # printf "oom found"
    gene=$(dx watch $job --get-stdout\
         --no-job-info -q --color off |\
         tail -n1 |\
        #  awk -v FS="(/|.ssv)" '{print $(NF-1)}') # ssv will catch from carrier files
         awk -v FS="(/|.csv)" '{print $(NF-1)}') # csv will catch from variant files
    printf "$job\t$gene\n"
  fi
done > my_oom_list.txt