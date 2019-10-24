#!/bin/bash

tar -czvf src.tar.gz *.c *.h makefile
scp src.tar.gz 21716194@ecm-ubl-006.uniwa.uwa.edu.au:~/
