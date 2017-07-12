#!/usr/bin/env sh

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 designfile responsefile"
    exit 1
fi
designs=$1
responses=$2
if ! [ -f "$1" ] ; then
   echo "designfile $designs does not exist"
   exit 1
fi
cat $designs | parallel --verbose batch.py {} > $responses.tmp
process.py $designs $responses.tmp > $responses
cat $responses | awk '{print $NF}'

