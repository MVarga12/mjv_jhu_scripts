" vim dotfile of Matt Varga

" Macros and functions
    function! RenameFile()
        let old_name = expand('%')
        let new_name = input('New file name: ', expand('%'), 'file')
        if new_name != '' && new_name != old_name
            exec ':saveas ' . new_name
            exec ':silent !rm ' . old_name
            redraw!
        endif
    endfunction
    map <leader>n :call RenameFile()<cr>

    " automatically reload vimrc upon buffer write
    if has ('autocmd') " Remain compatible with earlier versions
        augroup vimrc     " Source vim configuration upon save
            autocmd! BufWritePost $MYVIMRC source % | echom "Reloaded " . $MYVIMRC | redraw
            autocmd! BufWritePost $MYGVIMRC if has('gui_running') | so % | echom "Reloaded " . $MYGVIMRC | endif | redraw
        augroup END
    endif " has autocmd

" Key bindings (vim generic, not package specific)
    "This unsets the "last search pattern" register by hitting return
        nnoremap <CR> :noh<CR><CR>

    " cpp comment line, also works on multiple lines in visual mode
        map <C-C> :TComment<CR>
    
    " toggle folds with space
        nnoremap <space> za 

    " jump to end of line and add a semi-colon in insert mode
        inoremap <Leader>; <C-o>A;
    
    " quicksave with ,s
        nnoremap <Leader>s :w<CR>

    " easier switching from any split (specifically term split)
        nnoremap <C-h> <C-w>h
        nnoremap <C-j> <C-w>j
        nnoremap <C-k> <C-w>k
        nnoremap <C-l> <C-w>l

    " Abbreviations (Text expansion)
        iab ilist <TAB>\begin{itemize} <CR>\end{itemize}
        iab ieq <TAB>\begin{equation*} <CR>\end{equation*}
        iab iaq <TAB>\begin{align*} <CR>\end{align*}
        iab idq <TAB>\begin{displayquote} <CR>\end{displayquote}
        iab vv std::array<double, 3><SPACE>
 
" Random other settings
    " open with folds
        set foldmethod=indent
        set foldlevelstart=1 " start with most folds open
        set foldnestmax=10 " limit nested folds
        "space toggles folds
    
    " Make folder for swap files
        set swapfile
        set dir=~/tmp

    set updatetime=100
    set autoread
    "set lazyredraw
    set ttyfast

    "split-term
        set nosplitright
        set splitbelow
    
" GUI Specific Settings
    set conceallevel=2
    set concealcursor=nvc
    let g:tex_conceal="adgms"

    " set font to Adobe Source Code Pro
        set gfn=Hasklig\ Medium\:h13
    
    " set colour scheme
        "set t_Co=256
        if has('nvim')
            let $NVIM_TUI_ENABLE_TRUE_COLOR=1
            let base16colorspace=256
            set background=dark
            set termguicolors "nvim
            colorscheme Tomorrow-Night-Eighties 
        else
            colorscheme desert
        endif
    
    " show line numbers
        set nu
        set showmode
    
    " dummy sign to key sign column permanently open
        autocmd BufEnter * sign define dummy
        autocmd BufEnter * execute 'sign place 9999 line=1 name=dummy buffer=' . bufnr('')
        
    " show leader in the bottom right hand corner
        set showcmd
    
    " set tab to four spaces
        set tabstop=4
        set shiftwidth=4
        set softtabstop=4
        set expandtab
    
    " Comments as italics (make sure terminfo has sitm="\E[3m" and ritm="\E[23m"
        let t_ZH="\e[3m"
        let t_ZR="\e[23m"
        highlight Comment gui=italic cterm=italic term=italic

    " Some custom syntax highlighting
        hi TodoPriorityHigh cterm=italic ctermfg=52 gui=italic guifg=#DB0700
        hi TodoPriorityLow cterm=italic ctermfg=52 gui=italic guifg=#DBC200
        hi TodoPriorityMed cterm=italic ctermfg=52 gui=italic guifg=#DB8700
        hi TodoDone cterm=italic ctermfg=52 gui=italic guifg=#00ACDB
        hi TodoInProgress cterm=italic ctermfg=52 gui=italic guifg=#9D00D8
        hi CodeNote cterm=italic ctermfg=52 gui=italic guifg=#00D692
        call matchadd('TodoPriorityHigh', '\s\=[pP][rR][iI][oO][rR][iI][tT][yY]\s[hH][iI][gG][hH]')
        call matchadd('TodoPriorityLow', '\s\=[pP][rR][iI][oO][rR][iI][tT][yY]\s[lL][oO][wW]')
        call matchadd('TodoPriorityMed', '\s\=[pP][rR][iI][oO][rR][iI][tT][yY]\s[mM][eE][dD][iI]\=[uU]\=[mM]\=')
        call matchadd('TodoInProgress', '\s\=[iI][nN]\s[pP][rR][oO][gR][rR][eE][sS][sS]')
        call matchadd('TodoDone', '\s\=D[oO][nN][eE]')
        call matchadd('CodeNote', '\s\=N[oO][tT][eE]')

    " Status Line
    " modified from http://got-ravings.blogspot.com/2008/08/vim-pr0n-making-statuslines-that-own.html
        set statusline=
        set statusline+=%f
        set statusline+=%m
        set statusline+=%R%=
        set statusline+=%< " folding left
        set statusline+=%{synIDattr(synID(line('.'),col('.'),1),'name')}\ " highlight type on word
        set statusline+=%{gutentags#statusline('[',']')}\ 
        set statusline+=%(%3l,%02c%03V%)\ " row,column,virtual-column
        "set statusline+=\b\:%-04O\ " cursor hex offset from start of file
        "set statusline+=\c\:%03b\ " char byte code under cursor
        set statusline+=[%n%W\,%{strlen(&ft)?&ft:'none'}]\ " flags and filetype
        set statusline+=[%p%%] " percentage of the file
