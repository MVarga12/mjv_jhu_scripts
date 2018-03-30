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

" Plugin options
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

    " ctrlp
    let g:ctrlp_match_window = 'bottom,order:btt,min:1,max:30,results:30'

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

" new blank tab 
    nnoremap <silent> <leader>tb :tabedit<CR>

" dummy sign to key sign column permanently open
    autocmd BufEnter * sign define dummy
    autocmd BufEnter * execute 'sign place 9999 line=1 name=dummy buffer=' . bufnr('')

" show tab numbers in tabline
" Rename tabs to show tab number.
" (Based on http://stackoverflow.com/questions/5927952/whats-implementation-of-vims-default-tabline-function)
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
