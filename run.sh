#!/bin/bash

export JAVA_HOME=/home/work/software/java
export HBASE_HOME=/home/work/software/hbase
export HADOOP_HOME=/home/work/software/hadoop

export PATH=$JAVA_HOME/bin:$HBASE_HOME/bin:$HADOOP_HOME/bin:$PATH
export LANG="en_US.UTF-8"

task_name="graph"
input=/user/work/zsl/flow_feature/flow/*
#input=/user/work/zsl/$task_name/test.csv
job_name="zhangshunlong_"$task_name

hadoop fs -rmr /user/work/zsl/$task_name/res

hadoop jar /home/work/software/hadoop/contrib/streaming/hadoop-streaming.jar \
    -D mapred.map.tasks=100 \
    -D mapred.reduce.tasks=20 \
    -D mapred.job.name=$job_name \
    -D stream.num.map.output.key.fields=2 \
    -D mapred.text.key.partitioner.options=-k1,1 \
    -input  $input \
    -output /user/work/zsl/$task_name/res \
    -mapper chouqu_mapper.py \
    -reducer chouqu_reducer.py \
    -file /home/zsl/tasks/$task_name/chouqu_mapper.py \
    -file /home/zsl/tasks/$task_name/chouqu_reducer.py \
    -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner

rm -rf res
hadoop fs -cat /user/work/zsl/data_info/res/* | sort -n -t $'\t' -k1,1 * -o res.txt