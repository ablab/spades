#!/usr/bin/env csh
echo "starting command "
nohup $argv[1] $argv[2] $argv[3] &
ps -f -U mchaisso 
