" set runtine path
set runtimepath^=~/.vim runtimepath+=~/.vim/after
let &packpath = &runtimepath

" filetype stuff
filetype on
filetype plugin on
filetype plugin indent on
autocmd BufNewFile,BufReadPost *.md set ft=markdown

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
        Plug 'tpope/vim-markdown'
	    Plug 'octol/vim-cpp-enhanced-highlight'
	    Plug 'dag/vim-fish', {'for' : 'fish'}

	    " Deoplete
	    if has('nvim')
	    	Plug 'Shougo/deoplete.nvim', {'do': ':UpdateRemotePlugins'}
	    	Plug 'zchee/deoplete-clang',
	    	Plug 'zchee/deoplete-jedi',
	    endif

	    " Linters
	    Plug 'Shougo/neoinclude.vim'
	    Plug 'neomake/neomake'
	    Plug 'lervag/vimtex'

	    " Formatting
	    Plug 'tomtom/tcomment_vim' " universal commenter for embedded filetypes
	    Plug 'rhysd/vim-clang-format' " plugin for clang-format
	    " Plug 'Yggdroot/indentLine' " shows indent level
	    Plug 'Raimondi/delimitMate' " autoclosing of delimiters

	    " Debugging
	    Plug 'rizzatti/dash.vim'

	    " Distraction free writing
	    Plug 'junegunn/goyo.vim', {'for' : ['tex','text','md']}
	    Plug 'junegunn/limelight.vim', {'for' : ['tex','text','md']}

	    " Better search tools
	    Plug 'junegunn/fzf.vim'

	    " Gui plugins
	    Plug 'kshenoy/vim-signature'
	    Plug 'scrooloose/nerdtree'
	    Plug 'ap/vim-buftabline'

	    " Bookmarks and Tags
	    "Plug 'MattesGroeger/vim-bookmarks'
	    if ((&ft != "tex") && (&ft != "text") && (&ft != "markdown"))
	    	Plug 'ludovicchabant/vim-gutentags'
	    	Plug 'majutsushi/tagbar'
	    endif
	    Plug 'Valloric/ListToggle'

	    " Automatic session creation
	    Plug 'tpope/vim-obsession'
	    Plug 'dhruvasagar/vim-prosession'

	    " Easier movement
	    Plug 'chaoren/vim-wordmotion' " more intuitive defintions of words for w and e movement
	    Plug 'wesQ3/vim-windowswap' " swap splits while retaining layout
	    " Plug 'matze/vim-move'
	    " Plug 'easymotion/vim-easymotion'
	    "Plug 'terryma/vim-multiple-cursors' " allows SublimeText-like multiple cursors with regex support

	    " Tmux stuff
	    if exists("$TMUX")
	    	Plug 'christoomey/vim-tmux-navigator'
	    	Plug 'tmux-plugins/vim-tmux-focus-events'
	    endif

	    " Git Integration (off right now, using grv)
	    " Plug 'tpope/vim-fugitive'
	    " Plug 'tpope/vim-rhubarb'

	    " Other
	    Plug 'machakann/vim-sandwich' " more 'vim-like' version of tpope/vim-surround
	    Plug 'mbbill/undotree' " graphical visualization of the vim undotree
	    " Plug 'tpope/vim-surround'
	call plug#end()

	" Plugin options
	" vim-procession
	let g:procession_tmux_title = 1

	" vim-fish
	autocmd Filetype fish compiler fish
	autocmd Filetype fish setlocal foldmethod=expr

	" clang-format
	"map <C-I> :python /Users/mvarga/.vim/bundle/vim-clang-format/clang-format.py<cr>
	"imap <C-I> <c-o>:python /Users/mvarga/.vim/bundle/vim-clang-format/clang-format.py<cr>

	let g:clang_format#code_style = 'Google'
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
				\"ColumnLimit" : "100",
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

	" NERDTree
	nnoremap <silent> <F7> :NERDTreeToggle <Enter> :wincmd p <CR>
	" let NERDTreeAutoDeleteBuffer = 1
	let NERDTreeMinimalUI = 1
	let NERDTreeDirArrows = 1
	let NERDTreeQuitOnOpen = 1
	let NERDTreeShowLineNumbers = 0

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

	" CPP Enhanced Highlighting
	let g:class_decl_highlight = 1
	let g:cpp_member_variable_highlight = 1
	let g:cpp_experimental_template_highlight = 1
	let g:cpp_class_scope_highlight = 1

	" UndoTree
	if has("persistent_undo")
		set undodir=~/.undodir/
		set undofile
		set undolevels=1000
		set undoreload=1000
		" let g:undotree_WindowLayout = 2
		let g:undotree_SetFocusWhenToggle = 1
		nmap <silent> <F5> :UndotreeToggle<CR>
	endif

	" Goyo
	autocmd! User GoyoEnter Limelight
	autocmd! User GoyoLeave Limelight!
	nmap <F6> :Goyo<CR>

	" changing from the default 80 to accomodate for UndoTree panel
	" let g:goyo_width = 120

	" Gutentags
	let g:gutentags_ctags_extra_args = ['--extra=+p','--fields=+iaS']

	" Tagbar
	autocmd FileType c,cpp,py,h,hpp nested :call tagbar#autoopen(1)
	nmap <silent> <F8> :TagbarToggle<CR>
	let g:tagbar_width = 50
	let g:tagbar_autoshowtag = 1 " automatically open folds to show tag

	" deoplete
	let g:deoplete#sources#clang#libclang_path='/usr/local/Cellar/llvm/5.0.0/lib/libclang.dylib'
	let g:deoplete#sources#clang#clang_header='/usr/local/Cellar/llvm/5.0.0/include/c++'
	let g:deoplete#enable_at_startup=1
	let g:deoplete#max_list = 30
	call deoplete#custom#option('smart_case', v:true)
	autocmd InsertLeave,CompleteDone * if pumvisible() == 0 | pclose | endif

	"Tab Complete
	inoremap <expr><tab> pumvisible() ? "\<c-n>" : "\<tab>"

	"neomake
	" normal mode (after 1s; no delay when writing).
	call neomake#configure#automake('rw',2500)
	let neomake_verbose = 1
	let g:neomake_latex_lacheck_maker = {
				\ 'exe' : 'pdflatex'
				\}

	" maker stuff
	let g:neomake_cpp_enabled_makers = ['clang', 'clangtidy']

	if (!isdirectory('src'))
		let g:neomake_cpp_clang_maker = {
					\ 'exe': 'clang',
					\ 'args': ['-std=c++11','-Wall', '-Iinclude', '-I../include','-Wextra']
					\ }
	else
		let g:neomake_cpp_clang_maker = {
					\ 'exe': 'clang',
					\ 'args': ['-std=c++11','-Wall', '-Iinclude', '-I../include', '-Wextra']
					\ }
	endif

	let g:nemomake_cpp_clangtidy_maker = {
				\ 'exe' : 'clang-tidy',
				\ 'args' : ['-checks="-*","cppcoreguidelines-*","clang-analyzer-*"','--','-std=c++11','-Iinclude','-I../include']
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
	let g:neoinclude#max_processes = 20 " max number of include files processed
	" let g:neoinclude#ctags_commands = ''

	if exists('g:GtkGuiLoaded')
		call rpcnotify(1, 'Gui', 'Font', 'Envy Code R 13')
endif
source ~/.vimrc
