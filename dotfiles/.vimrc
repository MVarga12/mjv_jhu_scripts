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
    "execute pathogen#infect()
    "syntax on
    "filetype plugin indent on
    "call pathogen#helptags()
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

    "split-term
        set nosplitright
        set splitbelow
    
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

    " show tab numbers, buffer edited status, splits, and Obsession status in tabline
    " [Number][+ if edited and unsaved][file name for all splits, ',' delimited][Obsession Status ([$] for tracked, [S] for paused)]
    " off right now, trying to use buffers only
        "set showtabline=2
        "set tabline=%!MyTabLine()  " custom tab pages line
        function! MyTabLine()
            let s = ''
            " loop through each tab page
            for i in range(tabpagenr('$'))
                if i + 1 == tabpagenr()
                    let s .= '%#TabLineSel#'
                else
                    let s .= '%#TabLine#'
                endif

                if i + 1 == tabpagenr()
                    let s .= '%#TabLineSel#' " WildMenu
                else
                    let s .= '%#Title#'
                endif

                " set the tab page number (for mouse clicks)
                let s .= '%' . (i + 1) . 'T '
                " set page number string
                let s .= i + 1 . ''
                " get buffer names and statuses
                let n = ''  " temp str for buf names
                let m = 0   " &modified counter
                let buflist = tabpagebuflist(i + 1)

                " loop through each buffer in a tab
                for b in buflist
                    if getbufvar(b, "&buftype") == 'help'
                        " let n .= '[H]' . fnamemodify(bufname(b), ':t:s/.txt$//')
                    elseif getbufvar(b, "&buftype") == 'quickfix'
                        " let n .= '[Q]'
                    elseif getbufvar(b, "&modifiable")
                        let n .= fnamemodify(bufname(b), ':t') . ', ' " pathshorten(bufname(b))
                    endif

                    if getbufvar(b, "&modified")
                        let m += 1
                    endif
                endfor

                " let n .= fnamemodify(bufname(buflist[tabpagewinnr(i + 1) - 1]), ':t')
                let n = substitute(n, ', $', '', '')

                " add modified label
                if m > 0
                    let s .= '+'
                    " let s .= '[' . m . '+]'
                endif

                if i + 1 == tabpagenr()
                    let s .= ' %#TabLineSel#'
                else
                    let s .= ' %#TabLine#'
                endif

                " add buffer names
                if n == ''
                    let s.= '[New]'
                else
                    let s .= n
                endif

                let s .= '%{ObsessionStatus()}'

                " switch to no underlining and add final space
                let s .= ' '
            endfor

            let s .= '%#TabLineFill#%T'
                " right-aligned close button
                " if tabpagenr('$') > 1
                "   let s .= '%=%#TabLineFill#%999Xclose'
                " endif
            return s
        endfunction

" Plugin options
    " vim-procession
        let g:procession_tmux_title = 1

    " vim-fish
        autocmd Filetype fish compiler fish
        autocmd Filetype fish setlocal foldmethod=expr

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
        " show buffers in the tabline
        let g:airline#extensions#tabline#enabled = 1
        let g:airline#extensions#tabline#fnamemod = 1
        let g:airline#extensions#tagbar#enabled = 1

    " FZF
        " this is where fzf lives
        set rtp+=/usr/local/opt/fzf

        " ignore these when searching
        
        " use CtrlP keybindings for fzf
        map <C-p> :FZF 
        let g:fzf_action = {
            \ 'ctrl-t': 'tab split',
            \ 'ctrl-x': 'split',
            \ 'ctrl-v': 'vsplit'
        \ }

        " match fzf preview window colors to that of the vim colorscheme
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

        " let fzf remember things
        let g:fzf_history_dir = '~/.local/share/fzf_history'

        " opens in fullscreen with a preview window of where the string is in the file
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
        "set statusline+=%#warningmsg#
        "set statusline+=%{SyntasticStatuslineFlag()}
        "set statusline+=%*
        "
        "let g:syntastic_always_populate_loc_list = 1
        "let g:syntastic_auto_loc_list = 1
        "let g:syntastic_check_on_open = 1
        "let g:syntastic_check_on_wq = 0
        "
        "" Syntastic C++11 support
        "let g:syntastic_cpp_compiler = 'clang++'
        "let g:syntastic_cpp_compiler_options = ' -std=c++11 -stdlib=libc++'
        ""let g:syntastic_cpp_include_dirs=['include','../include']
    
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

        let g:undotree_WindowLayout = 2

        " using relative positioning instead
        "let g:undotree_CustomUndotreeCmd = 'vertical 32 new'
        "let g:undotree_CustomDiffpanelCmd= 'belowright 12 new'
        "let g:undotree_SetFocusWhenToggle = 1

    " Goyo
        autocmd! User GoyoEnter Limelight
        autocmd! User GoyoLeave Limelight!
        nmap <F6> :Goyo<CR>

    " Gutentags
        let g:gutentags_ctags_extra_args = ['--extra=+p','--fields=+iaS']

    " Tagbar
        autocmd VimEnter * nested :call tagbar#autoopen(1)
        autocmd FileType * nested :call tagbar#autoopen(0)
        autocmd BufEnter * nested :call tagbar#autoopen(0)
        nmap <F8> :TagbarToggle<CR>
        "let g:tagbar_autoclose = 1
        "let g:tagbar_autofocus = 1
        "let g:tagbar_compact = 1
        "let g:tagbar_previewwin_pos = "aboveleft"
        let g:tagbar_width = 50

    " vim-bookmarks
        let g:bookmark_auto_close = 1
        let g:bookmark_annotation_sign = '♪'

    " incsearch
        map / <Plug>(incsearch-forward)
        map ? <Plug>(incsearch-backward)
        map g/ <Plug>(incsearch-stay)
    
    " deoplete
        let g:deoplete#sources#clang#libclang_path='/usr/local/Cellar/llvm/5.0.0/lib/libclang.dylib'
        let g:deoplete#sources#clang#clang_header='/usr/local/Cellar/llvm/5.0.0/include/c++'
        let g:deoplete#enable_at_startup=1
        "let g:deoplete#complete_method = 'completefunc'
        let g:deoplete#max_list = 20
        let g:deoplete#enable_smart_case = 1
        "autoclose scratch
        autocmd InsertLeave,CompleteDone * if pumvisible() == 0 | pclose | endif
        "Tab Complete
        inoremap <expr><tab> pumvisible() ? "\<c-n>" : "\<tab>"
    
    "neomake
        "" normal mode (after 1s; no delay when writing).
        "call neomake#configure#automake('rnw',5000)
        "let neomake_verbose = 1
        ""let g:neomake_cpp_enabled_makers = ['clang']
        ""let g:neomake_cpp_clang_args = ['-std=c++11','-Wall', '-Iinclude', '-Wextra', '-Wno-sign-conversion']
        "let g:neomake_cpp_clang_maker = {
        "    \ 'exe': 'clang',
        "    \ 'args': ['-std=c++11','-Wall', '-Iinclude', '-I../include', '-Wextra']
        "\ }
        "let g:neomake_open_list = 0
        "let g:neomake_warning_sign = {
        "    \ 'text': 'W',
        "    \ 'texthl': 'WarningMsg',
        "\ }
        "let g:neomake_error_sign = {
        "    \ 'text': 'E',
        "    \ 'texthl': 'ErrorMsg',
        "\ }
        "let g:neomake_info_sign = {
        "    \ 'text': 'i',
        "    \ 'texthl': 'NeomakeInfoSign',
        "\}

    " neoinclude
    "    let g:neoinclude#max_processes = 10 " max number of include files processed
    "    let g:neoinclude#ctags_commands = '/usr/local/bin/ctags'
    
    " ale
    "    let g:ale_cpp_clang_options = '-std=c++11 -I include -Wall'
    "    let g:ale_cpp_gcc_options = '-std=c++11 -I include -Wall'
    "    let g:ale_lint_on_text_changed = 0
    "    let g:ale_lint_on_enter = 1
    "    let g:ale_sign_error = 'E'
    "    let g:ale_sign_warning = 'W'
    "    let g:ale_set_highlighs = 0
    
    "goyo
        " changing from the default 80 to accomodate for UndoTree panel
        let g:goyo_width = 120

