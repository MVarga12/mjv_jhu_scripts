## fish script
# Created: 29-Jun-2018 09:42:33 AM EDT
# Modified: 29-Jun-2018 09:54:44 AM EDT
# Created by: Matthew Varga
# Purpose: Tests if verification has been saved. If not, it creates a verification file.
##

function starthpc
    if test -e ~/.ssh/control:gateway2.marcc.jhu.edu:22:mvarga3@jhu.edu
        echo "Connecting with saved verification..."
        ssh gateway2.marcc.jhu.edu -l mvarga3@jhu.edu
    else
        ssh -XfNM gateway2.marcc.jhu.edu -l mvarga3@jhu.edu ; starthpc
    end
end
