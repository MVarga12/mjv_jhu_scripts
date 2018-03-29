" Macros
    " spaces around single operator, e.g. 1=1 -> 1 = 1
    "                   put cursor here -> ^
    let @d = 'i la ' 

    " spaces around double operator, e.g. 1==1 -> 1 == 1
    "                   put cursor here -> ^
    let @f = 'i lla '

    " cpp comment line, also works on multiple lines in visual mode
	map <C-C> I//<ESC>
	map <C-T> ^xx

" Abbreviations (Text expansion)
    iab ilist <TAB>\begin{itemize} <CR>\end{itemize}
    iab ieq <TAB>\begin{equation*} <CR>\end{equation*}
    iab iaq <TAB>\begin{align*} <CR>\end{align*}
    iab idq <TAB>\begin{displayquote} <CR>\end{displayquote}
    iab vv std::array<double, 3><SPACE>

" vim plugin handling with pathogen
" git clone url ~/.vim/bundle/PATH
    execute pathogen#infect()
    syntax on
    filetype plugin indent on

" lightline
    " relative path
    
    let g:lightline = {
        \ 'component_function': {
        \   'filename': 'LightLineFilename'
        \ }
        \ }
    function! LightLineFilename()
        return expand('%')
    endfunction
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
    set statusline+=%#warningmsg#
    set statusline+=%{SyntasticStatuslineFlag()}
    set statusline+=%*
    
    let g:syntastic_always_populate_loc_list = 1
    let g:syntastic_auto_loc_list = 1
    let g:syntastic_check_on_open = 1
    let g:syntastic_check_on_wq = 0
    
    " Syntastic C++11 support
    let g:syntastic_cpp_compiler = 'clang++'
    let g:syntastic_cpp_compiler_options = ' -std=c++11 -stdlib=libc++'
    "let g:syntastic_cpp_include_dirs=['include','../include']

"CPP Enhanced Highlighting
    "let g:cpp_member_variable_highlight = 1
    let g:cpp_class_scope_highlight = 1

" NERDTree
    nnoremap <leader>f :NERDTreeToggle<Enter>
    nnoremap <silent> <leader>v :NERTreeFind<Cr>
    let NERDTreeQuitOnOpen = 1

" YCM
    let g:ycm_global_ycm_extra_conf = "~/.vim/.ycm_extra_conf.py"
    set completeopt-=preview
"    let g:ycm_show_diagnostics_ui = 0

" show leader in the bottom right hand corner
    set showcmd

" set tab to four spaces
    set tabstop=4
    set shiftwidth=4
    set softtabstop=4
    set expandtab

" easier window navigation in split gvim/vim
    nnoremap <C-J> <C-W><C-J>
    nnoremap <c-K> <C-W><C-K>
    nnoremap <C-L> <C-W><C-L>
    nnoremap <C-H> <C-W><C-H>

" always open with folds by indent
    set foldmethod=indent
    
" map % to v% so when you press % it jumps to the closing bracket and selects
" all text between the two
"noremap % v%

" ----------- GUI SHIT -----------
" set font to Adobe Source Code Pro
    set gfn=Source\ Code\ Pro\ Light:h13

" set colour scheme
    "set t_Co=256
    "let $NVIM_TUI_ENABLE_TRUE_COLOR=1
    "set termguicolors "nvim
    "colorscheme badwolf
    "colorscheme 256_noir
    "colorscheme desert
    "colorscheme apprentice
    colorscheme dracula
    "colorscheme gruvbox
    "set background=dark
    " Solarized
    "    set background=dark
    "    colorscheme Neosolarized

"This unsets the "last search pattern" register by hitting return
    nnoremap <CR> :noh<CR><CR>

" show line numbers
    set nu

" tagbar
    nnoremap <silent> <leader>tb :TagbarToggle<CR>

" dummy sign to key sign column permanently open
    autocmd BufEnter * sign define dummy
    autocmd BufEnter * execute 'sign place 9999 line=1 name=dummy buffer=' . bufnr('')
