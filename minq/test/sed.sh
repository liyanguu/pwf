#!/bin/bash
#filename: sed.sh
#usage: remove none-C chars

sed 's/[a-zA-Z0-9_{}();,:\[\] \t\n\*\/+-\#<>=!]/!d' $*
