set PATH $HOME/.jenv/bin $PATH
function fish_mode_prompt; end # turns off the vi indicator
set -g fish_key_bindings fish_vi_key_bindings
set -gx term $TERM screen-256color

if not set -q abbrs_initialized
    set -U abbrs_initialized
    echo -n Setting abbreviations...

    #COMMAND LINE ALIASES
    abbr lt "ls -lt | less"
    
    #MARCC ALIASES
    abbr starthpc "ssh -X \"mvarga3@jhu.edu\"@gateway2.marcc.jhu.edu"
    abbr start_gateway "ssh -fNM gateway2.marcc.jhu.edu -l mvarga3@jhu.edu"
    abbr stop_gateway "ssh -O stop gateway2.marcc.jhu.edu -l mvarga3@jhu.edu"
    
    #PROGRAM ALIASES
    abbr rendervmd "/Applications/VMD\ 1.9.3.app/Contents/vmd/tachyon_MACOSXX86 -aasamples 12 vmdscene.dat -format TARGA -o vmdscene.dat.tga"
    abbr prev "open -a Preview"
    abbr fixbox "epstool --addtiff4-preview --bbox"
    abbr rulebender "/Applications/RuleBender-2.2.1-osx64/RuleBender.app/Contents/MacOS/RuleBender ; exit;"
    abbr vmd "/Applications/VMD\ 1.9.3.app/Contents/MacOS/startup.command ; exit;"
    abbr rsy "rsync -r -p -t -l -D -u -v"
    abbr cpr "cp -r"
    abbr back "cd -"
    abbr run_bng "/Applications/RuleBender-2.2.1-osx64/BioNetGen-2.3/BNG2.pl"
    abbr pathogen_dir "cd ~/.vim/bundle"
    
    echo 'Done'
end

function fish_greeting
    fortune -a
end

funcsave fish_greeting

test -e {$HOME}/.iterm2_shell_integration.fish ; and source {$HOME}/.iterm2_shell_integration.fish
