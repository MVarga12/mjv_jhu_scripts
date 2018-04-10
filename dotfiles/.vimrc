" vim dotfile of Matt Varga

" Macros
    " spaces around single operator, e.g. 1=1 -> 1 = 1
    "                   put cursor here -> ^
    let @d = 'i la ' 

    " spaces around double operator, e.g. 1==1 -> 1 == 1
    "                   put cursor here -> ^
    let @f = 'i lla '

" Key bindings (vim generic, not package specific)
    "This unsets the "last search pattern" register by hitting return
        nnoremap <CR> :noh<CR><CR>

    " cpp comment line, also works on multiple lines in visual mode
    	map <C-C> I//<ESC>
    	map <C-T> ^xx
    
     nnoremap <space> za 

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

" vim plugin handling with pathogen
    execute pathogen#infect()
    syntax on
    filetype plugin indent on
    call pathogen#helptags()
    "call pathogen#infect("after")
 
" Random other settings
    " open with folds
        set foldmethod=indent
        set foldlevelstart=10 " start with most folds open
        set foldnestmax=10 " limit nested folds
        "space toggles folds
    
    " Make folder for swap files
        set swapfile
        set dir=~/tmp

    set updatetime=100
    set autoread
    set lazyredraw
    set ttyfast
    
" GUI Specific Settings
    " set font to Adobe Source Code Pro
        set gfn=Hasklig\ Medium\:h13
    
    " set colour scheme
        "set t_Co=256
        let $NVIM_TUI_ENABLE_TRUE_COLOR=1
        let base16colorspace=256
        set background=dark
        set termguicolors "nvim
    
        " Pencil
        "    colorscheme pencil
        "    let g:pencil_terminal_italics = 1
        "    set background=dark
        
        "colorscheme deus
        colorscheme Tomorrow-Night-Eighties 
        "colorscheme badwolf
        "colorscheme hybrid
        " Solarized
        "    set background=dark
        "    colorscheme Neosolarized
    
    " show line numbers
        set nu
        set noshowmode
        "set cursorline
    
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

    "" show tab numbers in tabline
    "" Rename tabs to show tab number.
    "" (Based on http://stackoverflow.com/questions/5927952/whats-implementation-of-vims-default-tabline-function)
        if exists("+showtabline")
            function! MyTabLine()
                let s = ''
                let wn = ''
                let t = tabpagenr()
                let i = 1
                while i <= tabpagenr('$')
                    let buflist = tabpagebuflist(i)
                    let winnr = tabpagewinnr(i)
                    let s .= '%' . i . 'T'
                    let s .= (i == t ? '%1*' : '%2*')
                    let s .= ' '
                    let wn = tabpagewinnr(i,'$')
                    let s .= '%#TabNum#'
                    let s .= i
                    " let s .= '%*'
                    let s .= (i == t ? '%#TabLineSel#' : '%#TabLine#')
                    let bufnr = buflist[winnr - 1]
                    let file = bufname(bufnr)
                    let buftype = getbufvar(bufnr, 'buftype')
                    if buftype == 'nofile'
                        if file =~ '\/.'
                            let file = substitute(file, '.*\/\ze.', '', '')
                        endif
                    else
                        let file = fnamemodify(file, ':p:t')
                    endif
                    if file == ''
                        let file = '[No Name]'
                    endif
                    let s .= ' ' . file . ' '
                    let i = i + 1
                endwhile
                let s .= '%T%#TabLineFill#%='
                let s .= (tabpagenr('$') > 1 ? '%999XX' : 'X')
                return s
            endfunction
            set stal=2
            set tabline=%!MyTabLine()
            set showtabline=1
            highlight link TabNum Special
        endif

" Plugin options
    " clang-format
        "map <C-I> :python /Users/mvarga/.vim/bundle/vim-clang-format/clang-format.py<cr>
        "imap <C-I> <c-o>:python /Users/mvarga/.vim/bundle/vim-clang-format/clang-format.py<cr>
        
        let g:clang_format#code_style = 'Webkit'
        let g:clang_format#style_options = {
                    \"Standard" : "C++11",
                    \"BreakBeforeBraces" : "Custom",
                    \"BraceWrapping" : {
                    \    "BeforeCatch" : "true",
                    \    "AfterStruct" : "false",
                    \    "AfterClass" : "false",
                    \    "BeforeElse" : "true",
                    \    "SplitEmptyRecord" : "false",
                    \    "SplitEmptyFunction" : "false"
                    \    },
                    \"Language" : "Cpp",
                    \"AlwaysBreakTemplateDeclarations" : "true",
                    \"AlignOperands" : "true",
                    \"AlignTrailingComments" : "true",
                    \"AllowShortBlocksOnASingleLine" : "false",
                    \"AllowShortFunctionsOnASingleLine" : "false",
                    \"IncludeBlocks" : "Regroup",
                    \"PointerAlignment" : "Right",
                    \"ReflowComments" : "true",
                    \"ColumnLimit" : "200",
                    \"SpaceBeforeAssignmentOperators" : "true",
                    \"SpaceBeforeParens" : "ControlStatements",
                    \"Cpp11BracedListStyle" : "true",
                    \"SpacesInParentheses" : "false",
                    \"SpacesInSquareBrackets" : "false"
                    \}

        let g:clang_format#command = "clang-format"
        autocmd Filetype c,cpp ClangFormatAutoEnable

    " Airline
    let g:airline#extensions#tagbar#enabled = 1

    " FZF
        set rtp+=/usr/local/opt/fzf
        
        " use CtrlP keybindings
        map <C-p> :FZF 
        let g:fzf_action = {
            \ 'ctrl-t': 'tab split',
            \ 'ctrl-x': 'split',
            \ 'ctrl-v': 'vsplit'
        \ }

        let g:fzf_colors = {
            \'fg':       ['fg', 'Normal'],
            \ 'bg':      ['bg', 'Normal'],
            \ 'hl':      ['fg', 'Comment'],
            \ 'fg+':     ['fg', 'CursorLine', 'CursorColumn', 'Normal'],
            \ 'bg+':     ['bg', 'CursorLine', 'CursorColumn'],
            \ 'hl+':     ['fg', 'Statement'],
            \ 'info':    ['fg', 'PreProc'],
            \ 'border':  ['fg', 'Ignore'],
            \ 'prompt':  ['fg', 'Conditional'],
            \ 'pointer': ['fg', 'Exception'],
            \ 'marker':  ['fg', 'Keyword'],
            \ 'spinner': ['fg', 'Label'],
            \ 'header':  ['fg', 'Comment']
        \}

        let g:fzf_history_dir = '~/.local/share/fzf_history'
        command! -bang -nargs=* Ag
            \ call fzf#vim#ag(<q-args>,
            \                 <bang>0 ? fzf#vim#with_preview('up:60%')
            \                         : fzf#vim#with_preview('right:50%:hidden', '?'),
            \                 <bang>0)

        " Some key mappings
            noremap <Leader>b :Buffers<CR>
            noremap <Leader>h :History<CR>
            noremap <Leader>f :Ag<CR>

    " lightline
        "let g:lightline = {
        "    \ 'colorscheme': 'jellybeans',
        " \ }

        "" relative path
        "let g:lightline = {
        "    \ 'component_function': {
        "    \   'filename': 'LightLineFilename'
        "    \ }
        "    \ }
        "function! LightLineFilename()
        "    return expand('%')
        "endfunction

    " vimtex
        " make it so quickfix window does not open if only warnings, no errors
        let g:vimtex_quickfix_open_on_warning = 0
        let g:vimtex_quickfix_latexlog = {
            \ 'overfull' : 0,
            \ 'underful' : 0,
            \ 'packages' : {
                \ 'default' : 0,
            \},
        \}
    
    " Syntastic
    "    set statusline+=%#warningmsg#
    "    set statusline+=%{SyntasticStatuslineFlag()}
    "    set statusline+=%*
    "    
    "    let g:syntastic_always_populate_loc_list = 1
    "    let g:syntastic_auto_loc_list = 1
    "    let g:syntastic_check_on_open = 1
    "    let g:syntastic_check_on_wq = 0
    "    
    "    " Syntastic C++11 support
    "    let g:syntastic_cpp_compiler = 'clang++'
    "    let g:syntastic_cpp_compiler_options = ' -std=c++11 -stdlib=libc++'
    "    "let g:syntastic_cpp_include_dirs=['include','../include']
    
    " CPP Enhanced Highlighting
        let g:class_decl_highlight = 1
        let g:cpp_member_variable_highlight = 1
        let g:cpp_experimental_template_highlight = 1
        let g:cpp_class_scope_highlight = 1

    " UndoTree
        nmap <F5> :UndotreeToggle<CR>
        if has("persistent_undo")
            set undodir=~/.undodir/
            set undofile
        endif
        " using relative positioning instead
        let g:undotree_CustomUndotreeCmd = 'vertical 32 new'
        let g:undotree_CustomDiffpanelCmd= 'belowright 12 new'
        let g:undotree_SetFocusWhenToggle = 1

    " Goyo
        autocmd! User GoyoEnter Limelight
        autocmd! User GoyoLeave Limelight!
        nmap <F6> :Goyo<CR>

    " Tagbar
        nmap <F8> :TagbarToggle<CR>
        "let g:tagbar_autoclose = 1
        let g:tagbar_autofocus = 1
        let g:tagbar_compact = 1

    " vim-bookmarks
    let g:bookmark_highlight_lines = 1
    let g:bookmark_auto_close = 1
    let g:bookmark_annotation_sign = '♪'
