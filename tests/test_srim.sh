#!/usr/bin/bash

#     ./test_srim.sh <motion.zip>
# or
#      bash test_srim.sh <motion.zip>

python -m ssid srim $1 --ss-size 30 --arx-order 100 --inputs '[3, 6, 17]' --outputs '[12, 14, 19]'
