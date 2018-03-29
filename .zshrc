# If you come from bash you might have to change your $PATH.
# export PATH=$HOME/bin:/usr/local/bin:$PATH

# Tmuxinator

# Path to your oh-my-zsh installation.
export ZSH=/Users/mvarga/.oh-my-zsh

# Set name of the theme to load. Optionally, if you set this to "random"
# it'll load a random theme each time that oh-my-zsh is loaded.
# See https://github.com/robbyrussell/oh-my-zsh/wiki/Themes
ZSH_THEME=pygmalion
#ZSH_THEME="powerlevel9k/powerlevel9k"

# Uncomment the following line to use case-sensitive completion.
# CASE_SENSITIVE="true"

# Uncomment the following line to use hyphen-insensitive completion. Case
# sensitive completion must be off. _ and - will be interchangeable.
# HYPHEN_INSENSITIVE="true"

# Uncomment the following line to disable bi-weekly auto-update checks.
# DISABLE_AUTO_UPDATE="true"

# Uncomment the following line to change how often to auto-update (in days).
# export UPDATE_ZSH_DAYS=13

# Uncomment the following line to disable colors in ls.
# DISABLE_LS_COLORS="true"

# Uncomment the following line to disable auto-setting terminal title.
# DISABLE_AUTO_TITLE="true"

# Uncomment the following line to enable command auto-correction.
# ENABLE_CORRECTION="true"

# Uncomment the following line to display red dots whilst waiting for completion.
# COMPLETION_WAITING_DOTS="true"

# Uncomment the following line if you want to disable marking untracked files
# under VCS as dirty. This makes repository status check for large repositories
# much, much faster.
# DISABLE_UNTRACKED_FILES_DIRTY="true"

# Uncomment the following line if you want to change the command execution time
# stamp shown in the history command output.
# The optional three formats: "mm/dd/yyyy"|"dd.mm.yyyy"|"yyyy-mm-dd"
# HIST_STAMPS="mm/dd/yyyy"

# Would you like to use another custom folder than $ZSH/custom?
# ZSH_CUSTOM=/path/to/new-custom-folder

# Which plugins would you like to load? (plugins can be found in ~/.oh-my-zsh/plugins/*)
# Custom plugins may be added to ~/.oh-my-zsh/custom/plugins/
# Example format: plugins=(rails git textmate ruby lighthouse)
# Add wisely, as too many plugins slow down shell startup.
plugins=(git)

source $ZSH/oh-my-zsh.sh
source ~/.bin/tmuxinator.zsh

# User configuration

# export MANPATH="/usr/local/man:$MANPATH"

# You may need to manually set your language environment
# export LANG=en_US.UTF-8

# Preferred editor for local and remote sessions
# if [[ -n $SSH_CONNECTION ]]; then
#   export EDITOR='vim'
# else
#   export EDITOR='mvim'
# fi

# Compilation flags
# export ARCHFLAGS="-arch x86_64"

# ssh
# export SSH_KEY_PATH="~/.ssh/rsa_id"

# Set personal aliases, overriding those provided by oh-my-zsh libs,
# plugins, and themes. Aliases can be placed here, though oh-my-zsh
# users are encouraged to define aliases within the ZSH_CUSTOM folder.
# For a full list of active aliases, run `alias`.
#
# Example aliases
# alias zshconfig="mate ~/.zshrc"
# alias ohmyzsh="mate ~/.oh-my-zsh"

#COMMAND LINE ALIASES
alias lt="ls -lt"
alias lsdir="ls -p | grep \"/\" " 

#MARCC ALIASES
alias starthpc="ssh -X \"mvarga3@jhu.edu\"@gateway2.marcc.jhu.edu"
alias start_gateway="ssh -fNM gateway2.marcc.jhu.edu -l mvarga3@jhu.edu"
alias stop_gateway="ssh -O stop gateway2.marcc.jhu.edu -l mvarga3@jhu.edu"

#PROGRAM ALIASES
alias rendervmd="/Applications/VMD\ 1.9.3.app/Contents/vmd/tachyon_MACOSXX86 -aasamples 12 vmdscene.dat -format TARGA -o vmdscene.dat.tga"
alias prev="open -a Preview"
alias fixbox="epstool --addtiff4-preview --bbox"
alias rulebender="/Applications/RuleBender-2.2.1-osx64/RuleBender.app/Contents/MacOS/RuleBender ; exit;"
alias vmd="/Applications/VMD\ 1.9.3.app/Contents/MacOS/startup.command ; exit;"
alias rsy="rsync -r -p -t -l -D -u -v"
alias cpr="cp -r"
alias back="cd -"
alias run_bng="/Applications/RuleBender-2.2.1-osx64/BioNetGen-2.3/BNG2.pl"
alias pathogen_dir="cd ~/.vim/bundle"
alias gdat_plot="~/Google\ Drive/github_repos/jhu_scripts/nfsim_analysis_scripts/gdat_plot.py" #gdat_plot filename target
alias emacs="/usr/local/Cellar/emacs/25.3/Emacs.app/Contents/MacOS/Emacs"
alias emacsnw="/usr/local/Cellar/emacs/25.3/Emacs.app/Contents/MacOS/Emacs -nw"
