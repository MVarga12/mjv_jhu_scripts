" set runtine path
    set runtimepath^=~/.vim runtimepath+=~/.vim/after
    let &packpath = &runtimepath

    " filetype stuff
        filetype on
        filetype plugin on
        filetype plugin indent on

" map leader to comma
    let mapleader = ","

" read skeleton files for code headers
    if (&readonly == "off")
        let b:hfile="/Users/mvarga/.vim/headers/header.".expand("%:e")
        autocmd BufNewFile * let b:hfile="/Users/mvarga/.vim/headers/header.".expand("%:e") | call ReadHeader("hfile")
        autocmd BufRead * let b:hfile="/Users/mvarga/.vim/headers/header.".expand("%:e") | call AlterHeader("hfile")

        fun! ReadHeader(name)
            if filereadable(b:{a:name})
                augroup templates
                    au!
                    autocmd Filetype * silent! execute "0r $HOME/.vim/headers/header.".expand("%:e")
                    autocmd Filetype * %substitute#\[:VIM_EVAL_NEW:\]\(.\{-\}\)\[:END_EVAL_NEW:\]#\=eval(submatch(1))#ge
                    autocmd Filetype * %substitute#\[:VIM_EVAL_POST:\]\(.\{-\}\)\[:END_EVAL_POST:\]#\=eval(submatch(1))#ge
                    autocmd BufWritePre * ks|call Modified()|'s
                augroup make_exe
                    autocmd BufWritePre *.py,*.sh if !filereadable(expand('%')) | let b:is_new = 1 | endif
                    autocmd BufWritePost *.py,*.sh if get(b:, 'is_new', 0) | silent execute '!chmod +x %' | endif
                augroup END
            else 
                augroup templates
                    au!
                augroup END
            endif
        endfun!

        fun! AlterHeader(name)
            if filereadable(b:{a:name})
                augroup templates
                    au!
                    autocmd BufWritePre * ks | call Modified() | 's
                augroup END
            else
                augroup templates
                    au!
                augroup END
            endif
        endfun!
        
        fun! Modified()
            if line("$") > 20
                let l = 20
            else
                let l = line("$")
            endif
            exe "1," . l . "g/Modified: /s/Modified: .*/Modified: " . strftime('%a %e %b %Y %X %Z')
        endfun!
    endif

"vim-plug
call plug#begin('~/.vim/plugged')
    " Coding Plugins
        " Syntax enhancements
            Plug 'octol/vim-cpp-enhanced-highlight'
            Plug 'dag/vim-fish'

        " Snippets
            Plug 'SirVer/ultisnips'
            Plug 'honza/vim-snippets'

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
            Plug 'Shougo/neoinclude.vim'
            Plug 'neomake/neomake'
            Plug 'lervag/vimtex'

        " Formatting
            Plug 'tomtom/tcomment_vim' " universal commenter for embedded filetypes
            Plug 'rhysd/vim-clang-format' " plugin for clang-format
            Plug 'Yggdroot/indentLine' " shows indent level
            Plug 'Raimondi/delimitMate' " autoclosing of delimiters
    
        " Debugging
            Plug 'sakhnik/nvim-gdb' " integration of gdb within vim
            Plug 'rizzatti/dash.vim'

    " Distraction free writing
        "Plug 'junegunn/goyo.vim'
        "Plug 'junegunn/limelight.vim'

    " Better search tools
        Plug 'junegunn/fzf.vim'
        "Plug 'haya14busa/incsearch.vim'

    " Gui plugins
        Plug 'ap/vim-buftabline'

    " Bookmarks and Tags
        "Plug 'MattesGroeger/vim-bookmarks'
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
        Plug 'christoomey/vim-tmux-navigator'

    " Git Integration (off right now, using grv)
        Plug 'tpope/vim-fugitive'
        "Plug 'tpope/vim-rhubarb'

    " Other
        "Plug 'vim-scripts/timestamp.vim'
        Plug 'mbbill/undotree' " graphical visualization of the vim undotree

call plug#end()

" Plugin options
    " ultisnips
        let g:UltiSnipsExpandTrigger="<C-j>"
        let g:UltiSnipsJumpForwardTrigger="<C-b>"
        let g:UltiSnipsJumpBackwardTrigger="<C-z>"

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
                    \"AlignAfterOpenBracket" : "Align",
                    \"AlignOperands" : "true",
                    \"AlignTrailingComments" : "true",
                    \"AllowAllParametersOfDeclarationOnNextLine" : "true",
                    \"AllowShortBlocksOnASingleLine" : "false",
                    \"AllowShortIfStatementsOnASingleLine" : "false",
                    \"AllowShortLoopsOnASingleLine" : "false",
                    \"AlwaysBreakTemplateDeclarations" : "true",
                    \"BraceWrapping" : {
                    \    "AfterClass" : "false",
                    \    "AfterStruct" : "false",
                    \    "BeforeCatch" : "true",
                    \    "BeforeElse" : "true",
                    \    "SplitEmptyFunction" : "false",
                    \    "SplitEmptyRecord" : "false"
                    \    },
                    \"BreakBeforeBraces" : "Custom",
                    \"BreakConstructorInitializers" : "BeforeComma",
                    \"BreakStringLiterals" : "false",
                    \"ColumnLimit" : "200",
                    \"Cpp11BracedListStyle" : "true",
                    \"IncludeBlocks" : "Regroup",
                    \"Language" : "Cpp",
                    \"PenaltyBreakFirstLessLess" : "1000",
                    \"PointerAlignment" : "Right",
                    \"ReflowComments" : "true",
                    \"SpaceBeforeAssignmentOperators" : "true",
                    \"SpaceBeforeParens" : "ControlStatements",
                    \"SpacesInParentheses" : "false",
                    \"SpacesInSquareBrackets" : "false",
                    \"Standard" : "C++11"
                    \}

        let g:clang_format#command = "clang-format"
        autocmd Filetype c,cpp ClangFormatAutoEnable

    " Airline
        " show buffers in the tabline
        "let g:airline#extensions#tabline#enabled = 1
        "let g:airline#extensions#tabline#fnamemod = 1
        "let g:airline#extensions#tagbar#enabled = 1

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
            noremap <Leader>f :FZF<CR>
            noremap <Leader>F :Ag<CR>

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
        "let g:vimtex_quickfix_open_on_warning = 0
        "let g:vimtex_quickfix_latexlog = {
        "    \ 'overfull' : 0,
        "    \ 'underful' : 0,
        "    \ 'packages' : {
        "        \ 'default' : 0,
        "    \},
        "\}
    
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
        "autocmd! User GoyoEnter Limelight
        "autocmd! User GoyoLeave Limelight!
        "nmap <F6> :Goyo<CR>

        " changing from the default 80 to accomodate for UndoTree panel
            "let g:goyo_width = 120

    " Gutentags
        let g:gutentags_ctags_extra_args = ['--extra=+p','--fields=+iaS']

    " Tagbar
        autocmd VimEnter * nested :call tagbar#autoopen(1)
        "autocmd FileType * nested :call tagbar#autoopen(0)
        autocmd BufEnter * nested :call tagbar#autoopen(0)
        nmap <F8> :TagbarToggle<CR>
        "let g:tagbar_autoclose = 1
        "let g:tagbar_autofocus = 1
        "let g:tagbar_compact = 1
        "let g:tagbar_previewwin_pos = "aboveleft"
        let g:tagbar_width = 50

    " vim-bookmarks
        let g:bookmark_save_per_working_dir = 1
        let g:bookmark_auto_close = 1
        let g:bookmark_auto_save = 1
        let g:bookmark_annotation_sign = '♪'
        let g:bookmark_manage_per_buffer = 1

    " incsearch
        "map / <Plug>(incsearch-forward)
        "map ? <Plug>(incsearch-backward)
        "map g/ <Plug>(incsearch-stay)
    
    " deoplete
        let g:deoplete#sources#clang#libclang_path='/usr/local/Cellar/llvm/5.0.0/lib/libclang.dylib'
        let g:deoplete#sources#clang#clang_header='/usr/local/Cellar/llvm/5.0.0/include/c++'
        let g:deoplete#enable_at_startup=1
        " let g:deoplete#ignore_sources = ['tag']
        "let g:deoplete#complete_method = 'completefunc'
        let g:deoplete#max_list = 20
        call deoplete#custom#option('smart_case', v:true)
        "autoclose scratch
        autocmd InsertLeave,CompleteDone * if pumvisible() == 0 | pclose | endif
        "Tab Complete
        inoremap <expr><tab> pumvisible() ? "\<c-n>" : "\<tab>"
    
    "neomake
        " normal mode (after 1s; no delay when writing).
        call neomake#configure#automake('rnw',2500)
        let neomake_verbose = 1
        let g:neomake_latex_lacheck_maker = {
            \ 'exe' : 'pdflatex'
        \}
        let g:neomake_cpp_clang_maker = {
            \ 'exe': 'clang',
            \ 'args': ['-std=c++11','-Wall', '-Iinclude', '-I../include', '-Wextra']
        \ }
        let g:neomake_open_list = 0
        let g:neomake_warning_sign = {
            \ 'text': 'W',
            \ 'texthl': 'WarningMsg',
        \ }
        let g:neomake_error_sign = {
            \ 'text': 'E',
            \ 'texthl': 'ErrorMsg',
        \ }
        let g:neomake_info_sign = {
            \ 'text': 'i',
            \ 'texthl': 'NeomakeInfoSign',
        \}

        nmap <Leader><Leader>m :Neomake<CR>
        
    " neoinclude
        let g:neoinclude#max_processes = 10 " max number of include files processed
        "let g:neoinclude#ctags_commands = '/usr/local/bin/ctags'
    
    " ale
        "let g:ale_cpp_clang_options = '-std=c++11 -I include -Wall'
        "let g:ale_cpp_gcc_options = '-std=c++11 -I include -Wall'
        "let g:ale_lint_on_text_changed = 0
        "let g:ale_lint_on_enter = 1
        "let g:ale_sign_error = 'E'
        "let g:ale_sign_warning = 'W'
        "let g:ale_set_highlighs = 0
    
source ~/.vimrc
