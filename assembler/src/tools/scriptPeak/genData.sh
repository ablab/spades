#!/bin/bash
rm -rf data
java GetData -f distance_estimation < $1
