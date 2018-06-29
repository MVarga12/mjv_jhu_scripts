# function to change the filetype of all files of a particular filetype in
# a particular directory

function rename
    set __bash_args $argv
    for file in *.$__bash_args[1]
        mv -v -- "$file" (basename $file .$__bash_args[1]).$__bash_args[2]
    end
end
