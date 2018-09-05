set PATH $HOME/.jenv/bin $PATH
set PATH $HOME/.cargo/bin $PATH
set PKG_CONFIG_PATH $PKG_CONFIG_PATH /usr/local/lib/pkgconfig /usr/X11/lib/pkgconfig/ /opt/X11/lib/pkgconfig/
set PYTHONPATH $PYTHONPATH /usr/local/lib/python2.7/site-packages
set QTDIR $QTDIR /usr/local/opt/qt5

abbr startpluto "ssh -Y mvarga3@pluto.bph.jhu.edu"
abbr btx "byobu-tmux"
abbr see "nvim -R -u ~/.see_vimrc"
abbr lt "ls -alth | nvim -R -u ~/.see_vimrc -"
abbr rendervmd '/Applications/VMD\ 1.9.3.app/Contents/vmd/tachyon_MACOSXX86 -aasamples 12 vmdscene.dat -format TARGA -o vmdscene.dat.tga'
abbr prev "open -a Preview"

set -g fish_key_bindings fish_vi_key_bindings
#function fish_mode_prompt; end # turns off the vi indicator
#set -gx term $TERM screen-256color

# Less Colors for Man Pages
set -gx LESS_TERMCAP_mb \e'[01;31m'       # begin blinking
set -gx LESS_TERMCAP_md \e'[01;38;5;74m'  # begin bold
set -gx LESS_TERMCAP_me \e'[0m'           # end mode
set -gx LESS_TERMCAP_se \e'[0m'           # end standout-mode
set -gx LESS_TERMCAP_so \e'[38;5;246m'    # begin standout-mode - info box
set -gx LESS_TERMCAP_ue \e'[0m'           # end underline
set -gx LESS_TERMCAP_us \e'[04;38;5;146m' # begin underline

function fish_greeting
    fortune -a
end

funcsave fish_greeting

test -e {$HOME}/.iterm2_shell_integration.fish ; and source {$HOME}/.iterm2_shell_integration.fish
set -g fish_user_paths "/usr/local/opt/sqlite/bin" $fish_user_paths
set -g fish_user_paths "/usr/local/sbin" $fish_usr_paths
set -g fish_user_paths "/usr/local/opt/qt5/bin" $fish_user_paths
# set -g fish_user_paths "/usr/local/opt/qt@5.5/bin" $fish_user_paths
# set -g fish_user_paths "/usr/local/opt/qt5/bin" $fish_user_paths
