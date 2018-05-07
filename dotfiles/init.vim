" set runtine path
set runtimepath^=~/.vim runtimepath+=~/.vim/after
let &packpath = &runtimepath

" filetype stuff
filetype on
filetype plugin on
filetype plugin indent on
autocmd BufNewFile,BufReadPost *.md set ft=markdown
autocmd BufNewFile,BufReadPost *.mdp set ft=dosini " this is solely for commenting with t-comment

" map leader to comma
let mapleader = "\<SPACE>"

" read skeleton files for code headers
if (&readonly == "off")
    let b:hfile="/Users/mvarga/.vim/headers/header.".expand("%:e")
    autocmd BufNewFile * let b:hfile="/Users/mvarga/.vim/headers/header.".expand("%:e") | call ReadHeader("hfile")
    autocmd BufRead,BufEnter * let b:hfile="/Users/mvarga/.vim/headers/header.".expand("%:e") | call AlterHeader("hfile")

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
            exe "1," . l . "g/Modified: /s/Modified: .*/Modified: " . strftime('%v %I:%M:%S %p %Z')
        endfun!
endif

"vim-plug
    call plug#begin('~/.vim/plugged')
    " Coding Plugins
    " Syntax enhancements
    Plug 'tpope/vim-markdown', {'for' : ['mkd','md','markdown']}
    Plug 'octol/vim-cpp-enhanced-highlight'
    Plug 'dag/vim-fish', {'for' : 'fish'}

    " Linters
    " Plug 'Shougo/neoinclude.vim'
    " Plug 'neomake/neomake'
    " Plug 'lervag/vimtex'
    Plug 'w0rp/ale'
    Plug 'roxma/nvim-completion-manager'
    Plug 'roxma/ncm-clang'
    inoremap <expr> <Tab> pumvisible() ? "\<C-n>" : "\<Tab>"
    inoremap <expr> <S-Tab> pumvisible() ? "\<C-p>" : "\<S-Tab>"
    Plug 'autozimu/LanguageClient-neovim', {
                \ 'branch' : 'next',
                \ 'do' : 'bash install.sh'
                \ }

    " Formatting
    Plug 'tomtom/tcomment_vim' " universal commenter for embedded filetypes
    Plug 'rhysd/vim-clang-format' " plugin for clang-format

    " Debugging
    Plug 'rizzatti/dash.vim'

    " writing
    Plug 'reedes/vim-pencil'
    Plug 'reedes/vim-litecorrect'
    Plug 'junegunn/goyo.vim'
    Plug 'junegunn/limelight.vim'
    Plug 'chrisbra/unicode.vim'

    " Snippets
    Plug 'SirVer/ultisnips'
    Plug 'honza/vim-snippets'

    " Better search tools
    Plug 'junegunn/fzf.vim'


    " Bookmarks and Tags
    Plug 'MattesGroeger/vim-bookmarks'
    Plug 'ludovicchabant/vim-gutentags'
    Plug 'majutsushi/tagbar'
    Plug 'Valloric/ListToggle'

    " Automatic session creation
    Plug 'tpope/vim-obsession'
    Plug 'dhruvasagar/vim-prosession'

    " Easier movement
    Plug 'chaoren/vim-wordmotion' " more intuitive defintions of words for w and e movement
    Plug 'wesQ3/vim-windowswap' " swap splits while retaining layout

    " Tmux stuff
    if exists("$TMUX")
        Plug 'christoomey/vim-tmux-navigator'
        Plug 'tmux-plugins/vim-tmux-focus-events'
    endif

    " Other
    Plug 'mbbill/undotree' " graphical visualization of the vim undotree
call plug#end()

" Plugin options
    " vim-procession
    let g:procession_tmux_title = 1

    " vim-fish
    autocmd Filetype fish compiler fish
    autocmd Filetype fish setlocal foldmethod=expr

    " clang-format
    let col_width = 120 " column width for ClangFormat
    let g:clang_format#code_style = 'Webkit'
    let g:clang_format#style_options = {
                \"ColumnLimit" : col_width,
                \"AlwaysBreakTemplateDeclarations" : "true"
    \}

    let g:clang_format#command = "clang-format"
    autocmd Filetype c,cpp ClangFormatAutoEnable

    " FZF
        " this is where fzf lives
        set rtp+=/usr/local/opt/fzf

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

    " CPP Enhanced Highlighting
    let g:class_decl_highlight = 1
    let g:cpp_member_variable_highlight = 1
    let g:cpp_class_scope_highlight = 1

    " UndoTree
    if has("persistent_undo")
        set undodir=~/.undodir/
        set undofile
        set undolevels=1000
        set undoreload=1000
        let g:undotree_SetFocusWhenToggle = 1
        nmap <silent> <F5> :UndotreeToggle<CR>
    endif

    " Goyo
        autocmd! User GoyoEnter Limelight
        autocmd! User GoyoLeave Limelight!
        nmap <F6> :Goyo<CR>

        " changing from the default 80 to accomodate for UndoTree panel
        let g:goyo_width = 120

    " Gutentags
    let g:gutentags_ctags_extra_args = ['--extra=+p','--fields=+iaS']

    " Tagbar
    autocmd FileType c,cpp,py,h,hpp nested :call tagbar#autoopen(1)
    nmap <silent> <F8> :TagbarToggle<CR>
    let g:tagbar_width = 30
    let g:tagbar_autoshowtag = 1 " automatically open folds to show tag

    " Snippets
    let g:UltiSnipsExpandTrigger="<Tab>"
    let g:UltiSnipsJumpForwardTrigger="<Tab"
    let g:UltiSnipsJumpBackwardTrigger="<S-Tab>"

    " ALE
        " Generic Options
            " not using pylint for python, it's absolutely crazy
            let g:ale_linters = {'cpp': ['clang','clangtidy'], 'python': ['flake8']}
            let g:ale_lint_on_text_changed = 1
            let g:ale_lint_on_enter = 1
            let g:ale_sign_error = 'E'
            let g:ale_sign_warning = 'W'
            let g:ale_set_highlighs = 0
        "C++
            let g:ale_cpp_clang_options = '-std=c++11 -Iinclude -Wall'
            let g:ale_cpp_clangtidy_options = '-Wall -std=c++11 -x c++'
            let g:ale_c_build_dir = '.'
        " Python

    "some writing stuff
    " autocmd Filetype tex,latex,text,td,md,markdown,mkd
    autocmd BufAdd,BufReadPre *.tex,*.latex,*.tex,*.text,*.txt,*.md,*.markdown,*.mkd
                \ call pencil#init({'wrap':'soft','textwidth':120, 'conceallevel':2, 'autoformat':1})
                \ | call litecorrect#init()
                \ | setlocal spell spelllang=en_us noruler nonumber norelativenumber
                \ | setlocal foldopen+=search
                \ | setlocal nocursorcolumn

source ~/.vimrc
