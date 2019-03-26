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
        " nnoremap <Leader><Leader>s :source ~/.vimrc<CR>
        " if has ('autocmd') " Remain compatible with earlier versions
        "     augroup vimrc  " Source vim configuration upon save
        "         au!
        "         autocmd! BufWritePost $MYVIMRC source % | echom "Reloaded " . $MYVIMRC | redraw
        "         autocmd! BufWritePost ~/.vimrc if has('gui_running') | so % | echom "Reloaded " . $MYGVIMRC | endif | redraw
        "     augroup END
        " endif " has autocmd
    
    " make a scratch buffer which leaves no trace when it doesn't exist in the window anymore.
    " source: http://dhruvasagar.com/2014/03/11/creating-custom-scratch-buffers-in-vim
    function! ScratchEdit(cmd,options)
        exe a:cmd tempname()
        setlocal buftype=nofile bufhidden=wipe nobuflisted
        if !empty(a:options) | exe 'setlocal' a:options | endif
    endfunction!
    
    command! -bar -nargs=* Sedit call ScratchEdit('edit', <q-args>)
    command! -bar -nargs=* Ssplit call ScratchEdit('aboveleft 10split', <q-args>)
    command! -bar -nargs=* Svsplit call ScratchEdit('10vsplit', <q-args>)
    nnoremap <Leader>ss :Ssplit<CR>
    nnoremap <Leader>sv :Svsplit<CR>

" Key bindings (vim generic, not package specific)
    " capitalization key bindings
    " (http://vim.wikia.com/wiki/Capitalize_words_and_regions_easily)
        if (&tildeop)
            "capitalize word (cursor position to end of word)
            nmap gcw guw~l
            " capitalize WORD (cursor position to end of WORD)
            nmap gcW guW~l
            " capitalize inner word (start to end)
            nmap gciw guiw~l
            " capitalize inner WORD (start to end)
            nmap gciW guiW~l
            " capitalize inner sentence
            nmap gcis guis~l
            " capitalize until end of line (cursor position)
            nmap gc$ gu$~l
            " capitalize whole line (start to end)
            nmap gcgc guu~l
            " capitalize whole line
            nmap gcc guu~l
            " capitalize highlighted text
            vmap gc gu~l
        else
            nmap gcw guw~h
            nmap gcW guW~h
            nmap gciw guiw~h
            nmap gciW guiW~h
            nmap gcis guis~h
            nmap gc$ gu$~h
            nmap gcgc guu~h
            nmap gcc guu~h
            vmap gc gu~h
        endif
    "This unsets the "last search pattern" register by hitting return
        nnoremap <CR> :noh<CR><CR>
    
    " cpp comment line, also works on multiple lines in visual mode
        map <C-C> :TComment<CR>
    
    " jump to end of line and add a semi-colon in insert mode
        inoremap <Leader>; <C-o>A;
    
    " quicksave with ,s
        nnoremap <Leader>s :w<CR>

    " quickreload vimrc
        " nnoremap <Leader>r :source $MYVIMRC <CR>
    
    " save with sudo when opened without it
        cnoremap w!! w !sudo tee % > /dev/null
    
    " Easier resizing of panes
        nnoremap <Right> :vertical resize +2<CR>
        nnoremap <Left> :vertical resize -2<CR>
        nnoremap <Down> :resize +2<CR>
        nnoremap <Up> :resize -2<CR>
    
    " easier switching from any split (specifically term split)
        nnoremap <C-h> <C-w>h
        nnoremap <C-j> <C-w>j
        nnoremap <C-k> <C-w>k
        nnoremap <C-l> <C-w>l

    " make splits stay the same size when exiting one
        set noequalalways
    
    " Abbreviations (Text expansion)
        iab ilist <TAB>\begin{itemize} <CR>\end{itemize}
        iab ieq <TAB>\begin{equation*} <CR>\end{equation*}
        iab iaq <TAB>\begin{align*} <CR>\end{align*}
        iab idq <TAB>\begin{displayquote} <CR>\end{displayquote}
        iab vv std::array<double, 3><SPACE>

" Random other settings
    " open with folds
        set foldmethod=syntax

    " Make folder for swap files
        set swapfile
        set dir=~/tmp

        set updatetime=100
        set autoread
        set lazyredraw

    "split-term
        set nosplitright
        set splitbelow

    " spell check
        set spell spelllang=en_us

" GUI Specific Settings
    " press <Leader>z to zoom in on vim pane
        if exists('$TMUX')
            nnoremap <silent> <Leader>z :call system("tmux resize-pane -Z")<CR>
            " autocmd FocusGained * call system("tmux resize-pane -Z") " TODO: make it so that this doesn't fire if switching desktops
        endif

    " Focused pane numbering options 
        set number norelativenumber
        " augroup active_relative_number
        "     autocmd!
        "     autocmd BufEnter,FocusGained,InsertLeave * :setlocal number relativenumber
        "     autocmd BufLeave,FocusLost,InsertEnter * :setlocal number norelativenumber
        " augroup END

    " set colour scheme
    "set t_Co=256
        if has('nvim')
            let $NVIM_TUI_ENABLE_TRUE_COLOR=1
            let base16colorspace=256
            set termguicolors "nvim

            " TNE-Edited
            " set background=dark
            colorscheme pencil

           "Pencil Light
            " set background=light
            " colorscheme pencil
        else
            colorscheme desert
        endif
    
    " show line numbers
        set nu
        set relativenumber
        set showmode
    
    " dummy sign to key sign column permanently open
        autocmd BufEnter * sign define dummy
        autocmd BufEnter * execute 'sign place 9999 line=1 name=dummy buffer=' . bufnr('')

    " show leader in the bottom right hand corner
        set showcmd

    " Comments as italics (make sure terminfo has sitm="\E[3m" and ritm="\E[23m"
    " this is specifically for TMUX, but I don't want to wrap it in an if statement

    " Some custom syntax highlighting
        set list
        set listchars:trail:+
        " todos, their priorities, and notes syntax highlighting 
        function! SetItalicComments()
            let t_ZH="\e[3m"
            let t_ZR="\e[23m"
            if (&bg == "dark")
                highlight Comment gui=italic cterm=italic term=italic
            else 
                highlight Comment gui=italic cterm=italic term=italic
            endif
        endfunction!
        function! PriorityHighlighting()
            hi TodoPriorityHigh cterm=italic ctermfg=52 gui=italic guifg=#DB0700
            hi TodoPriorityLow cterm=italic ctermfg=52 gui=italic guifg=#DBC200
            hi TodoPriorityMed cterm=italic ctermfg=52 gui=italic guifg=#DB8700
            hi TodoDone cterm=italic ctermfg=52 gui=italic guifg=#00ACDB
            hi TodoInProgress cterm=italic ctermfg=52 gui=italic guifg=#9D00D8
            hi CodeNote cterm=italic ctermfg=52 gui=italic guifg=#00D692
            hi CodeIdea cterm=italic ctermfg=52 gui=italic guifg=#AF85FF
            hi LoneDash cterm=italic ctermfg=52 gui=italic guifg=#C3B3FF
            call matchadd('TodoPriorityHigh', '-\=\s\=[pP][rR][iI][oO][rR][iI][tT][yY]\s[hH][iI][gG][hH]:\=')
            call matchadd('TodoPriorityLow', '-\=\s\=[pP][rR][iI][oO][rR][iI][tT][yY]\s[lL][oO][wW]:\=')
            call matchadd('TodoPriorityMed', '-\=\s\=[pP][rR][iI][oO][rR][iI][tT][yY]\s[mM][eE][dD][iI]\=[uU]\=[mM]\=:\=')
            call matchadd('TodoInProgress', '\s\=[iI][nN]\s[pP][rR][oO][gR][rR][eE][sS][sS]:\=')
            call matchadd('TodoDone', '\s\=D[oO][nN][eE]:\=')
            call matchadd('CodeNote', '-*\s\=N[oO][tT][eE]:\=')
            call matchadd('CodeIdea', '-*\s\=I[dD][eE][aA]:\=')
            call matchadd('LoneDash', '\s-\s\=')
        endfunction!

        " ignore latex math in markdown
        function! MarkdownSyntaxIgnore()
            syntax region math start=/\$\$/ end=/\$\$/
            syntax match math '\$[^$].\{-}\$' 
            highlight link math Statement
        endfunction!

        highlight ColorColumn guibg = magenta 
        augroup highlight
            au!
            autocmd BufEnter,BufRead,BufNew * call SetItalicComments()
            autocmd BufEnter,BufRead,BufNew * call PriorityHighlighting()
            autocmd BufEnter,BufRead,BufNew *.md,*.mkd,*.markdown call MarkdownSyntaxIgnore()
            autocmd BufEnter, BufRead, BufNew * call matchadd('ColorColumn', '\%120v\S', 100)
        augroup END

    function! LinterStatus() abort
        let l:counts = ale#statusline#Count(bufnr(''))

        let l:all_errors = l:counts.error + l:counts.style_error
        let l:all_non_errors = l:counts.total - l:all_errors

        return l:counts.total == 0 ? 'OK' : printf(
        \   '%dW %dE',
        \   all_non_errors,
        \   all_errors
        \)
    endfunction

    " Status Line
    " modified from http://got-ravings.blogspot.com/2008/08/vim-pr0n-making-statuslines-that-own.html
    set statusline=
    set statusline+=%t
    set statusline+=%m
    set statusline+=%R%=
    set statusline+=%< " folding left
    set statusline+=%{synIDattr(synID(line('.'),col('.'),1),'name')}\ " highlight type on word
    set statusline+=%{gutentags#statusline('[',']')}\ 
    set statusline+=%(%3l,%02c%03V%)\ " row,column,virtual-column
    "set statusline+=\b\:%-04O\ " cursor hex offset from start of file
    "set statusline+=\c\:%03b\ " char byte code under cursor
    set statusline+=[\ALE:\ %{LinterStatus()}]\ 
    set statusline+=[%n%W\,%{strlen(&ft)?&ft:'none'}]\ " flags and filetype
    set statusline+=[%p%%] " percentage of the file

    " set tab to spaces
        set tabstop=8
        set softtabstop=4
        set shiftwidth=4
        set expandtab
        set smarttab

        " set tab to two spaces for markdown
        autocmd BufRead,BufEnter,BufNew *.md,*.mkd,*.markdown :setl shiftwidth=2
