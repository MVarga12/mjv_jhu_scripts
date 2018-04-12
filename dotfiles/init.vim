set runtimepath^=~/.vim runtimepath+=~/.vim/after
let &packpath = &runtimepath

"vim-plug
call plug#begin('~/.vim/plugged')
    " Coding Plugins
        " Syntax enhancements
            Plug 'octol/vim-cpp-enhanced-highlight'
            Plug 'dag/vim-fish'

        " Deoplete
            if has('nvim')
                Plug 'Shougo/deoplete.nvim', {'do': ':UpdateRemotePlugins'}
                Plug 'zchee/deoplete-clang' " Clang, C++, C completion library
                Plug 'zchee/deoplete-jedi'  " Python completion library
            else
                Plug 'Shougo/deoplete.nvim'
                Plug 'roxma/nvim-yarp'
                Plug 'roxma/vim-hug-neovim-rpc'
                Plug 'zchee/deoplete-clang'
            endif

        " Linters
            "Plug 'Shougo/neoinclude.vim'
            "Plug 'neomake/neomake'
            "Plug 'w0rp/ale'
            Plug 'lervag/vimtex'

        " Formatting
            Plug 'rhysd/vim-clang-format' " plugin for clang-format
            Plug 'Yggdroot/indentLine' " shows indent level
    
        " Debugging
            Plug 'sakhnik/nvim-gdb' " integration of gdb within vim

    " Distraction free writing
        Plug 'junegunn/goyo.vim'
        Plug 'junegunn/limelight.vim'

    " Better search tools
        "Plug 'usr/local/opt/fzf'
        Plug 'junegunn/fzf.vim'
        Plug 'haya14busa/incsearch.vim'

    " Gui plugins
        Plug 'vim-airline/vim-airline'
        Plug 'vim-airline/vim-airline-themes'

    " Bookmarks and Tags
        Plug 'MattesGroeger/vim-bookmarks'
        Plug 'ludovicchabant/vim-gutentags'
        Plug 'majutsushi/tagbar'
        Plug 'Valloric/ListToggle'

    " Automatic session creation
        Plug 'tpope/vim-obsession'
        Plug 'dhruvasagar/vim-prosession'

    " Easier movement
        Plug 'easymotion/vim-easymotion'
        Plug 'terryma/vim-multiple-cursors' " allows SublimeText-like multiple cursors with regex support
        Plug 'wesQ3/vim-windowswap' " swap splits while retaining layout

    " Git Integration (off right now, using grv)
        "Plug 'tpope/vim-fugitive'
        "Plug 'tpope/vim-rhubarb'

    " Other
        Plug 'vim-scripts/timestamp.vim'
        Plug 'mbbill/undotree' " graphical visualization of the vim undotree
call plug#end()

" All plugin specific commands are in vimrc
source ~/.vimrc
