#!/bin/bash

# testing each executable using input/hp1

$1/ad_sa1 ../../input/hp1 2 123456 | grep Final | awk -v p="ad_sa1 p1 2" -v t=1 -f utest.awk
$1/ad_sa2 ../../input/hp1 2 123456 | grep Final | awk -v p="ad_sa2 p1 2" -v t=1 -f utest.awk
$1/ad_rsa ../../input/hp1 2 123456 | grep Final | awk -v p="ad_rsa p1 2" -v t=1 -f utest.awk

# testing each executable using input/hp9

$1/ad_sa1 ../../input/hp9 2 123456 | grep Final | awk -v p="ad_sa1 p9 2" -v t=40 -f utest.awk
$1/ad_sa1 ../../input/hp9 3 123456 | grep Final | awk -v p="ad_sa1 p9 3" -v t=56 -f utest.awk
$1/ad_sa2 ../../input/hp9 2 123456 | grep Final | awk -v p="ad_sa2 p9 2" -v t=70 -f utest.awk
$1/ad_sa2 ../../input/hp9 3 123456 | grep Final | awk -v p="ad_sa2 p9 3" -v t=119 -f utest.awk
$1/ad_rsa ../../input/hp9 2 123456 | grep Final | awk -v p="ad_rsa p9 2" -v t=35 -f utest.awk

# EOF


