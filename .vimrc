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

" vim plugin handling with pathogen
" git clone url ~/.vim/bundle/PATH
    execute pathogen#infect()
    syntax on
    filetype plugin indent on

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
    nnoremap <C-K> <C-W><C-K>
    nnoremap <C-L> <C-W><C-L>
    nnoremap <C-H> <C-W><C-H>

" map % to v% so when you press % it jumps to the closing bracket and selects
" all text between the two
"noremap % v%

" ----------- GUI SHIT -----------
" set font to Adobe Source Code Pro
    set gfn=Source\ Code\ Pro\ Light:h13

" set colour scheme
    "colorscheme desert
    " Solarized
        "syntax enable "vim
        set termguicolors "nvim
        set background=dark
        colorscheme Neosolarized

"This unsets the "last search pattern" register by hitting return
    nnoremap <CR> :noh<CR><CR>

" show line numbers
    set nu

" tagbar
    nnoremap <silent> <leader>tb :TagbarToggle<CR>

" dummy sign to key sign column permanently open
    autocmd BufEnter * sign define dummy
    autocmd BufEnter * execute 'sign place 9999 line=1 name=dummy buffer=' . bufnr('')
