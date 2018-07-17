#!/usr/bin/env bash
##
# Created:  2-Jul-2018 10:21:58 AM EDT
# Modified:  2-Jul-2018 10:54:40 AM EDT
# Created by: Matthew Varga
# Purpose:
##

echo "Creating environment"
virtualenv ~/md_env
source ~/md_env/bin/activate

# Setup new python path
export PYTHONPATH_BK=$PYTHONPATH
export PYTHONPATH=~/md_env/lib/python2.7/site-packages

if [ "$PYTHONPATH" != "~/md_env/lib/python2.7/site-packages" ]; then
    echo "Exiting, python path not set up properly..."
    exit
fi

# install packages
echo "Installing required packages"
pip install --requirements md_env_reqs.txt
