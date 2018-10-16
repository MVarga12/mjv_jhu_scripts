# vmd tcl procedure: reload the last trajectory data set 
#                    for the current top molecule. 
#
# Copyright (c) 2004 by <Axel.Kohlmeyer@theochem.ruhr-uni-bochum.de>

proc reload {args} {
    set tmol [molinfo top]
    array set viewpoints {}
    foreach mol [molinfo list] {
    # save orientation and zoom parameters
      set viewpoints($mol) [molinfo $mol get { 
        center_matrix rotate_matrix scale_matrix global_matrix}]
    } 
    # delete all coordinates and (re)load the latest data set.
    animate delete all $tmol
    set files [lindex [molinfo $tmol get filename] 0]
    set lf [expr [llength $files] - 1]

    mol addfile [lindex $files $lf] mol $tmol \
        type [lindex [lindex [molinfo $tmol get filetype] 0] $lf]

    foreach mol [molinfo list] {
    # restore orientation and zoom
      molinfo $mol set {center_matrix rotate_matrix scale_matrix
        global_matrix} $viewpoints($mol)
    }

}
