#!/bin/bash

# Make a text file for plotting station names in GMT:
# x  y  size  angle  fontno  justify  text

# Get the lon, lat, and "text" from a station file:

#awk '{print $2, $3, "12", "0", "22", "CT", $1}'

awk '{print $2, $3, $1}'