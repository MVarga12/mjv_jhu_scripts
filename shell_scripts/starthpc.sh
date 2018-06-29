#!/usr/bin/env bash
##
# Created: 29-Jun-2018 09:55:17 AM EDT
# Modified: 29-Jun-2018 10:00:36 AM EDT
# Created by: Matthew Varga
# Purpose:
##

if [ -e "$HOME/.ssh/control:gateway2.marcc.jhu.edu:22:mvarga3@jhu.edu" ]
then
    echo "Starting MARCC with saved credentials..."
    ssh gateway2.marcc.jhu.edu -l mvarga3@jhu.edu
else
    echo "Creating gateway..."
    ssh -XfNM gateway2.marcc.jhu.edu -l mvarga3@jhu.edu
    ssh gateway2.marcc.jhu.edu -l mvarga3@jhu.edu
fi
