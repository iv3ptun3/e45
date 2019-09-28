#!/bin/bash

SIG=SIGINT

bjobs | \
    while read n u s q o; do
	if [ "$u" = "hayashu" ]; then
	    bkill -s $SIG $n &
	    usleep 10000
	fi
    done
