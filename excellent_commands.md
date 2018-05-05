# (n)vim commands
Note: all commands are in Normal mode or Command mode (denoted by ':'), unless otherwise stated
## Whitespace and Indenting
* `dw` - delete all whitespace online until a letter is reached
* `:%s/^\s\+//e` - deletes all leading whitespace in file
  * I think an easier way of doing this is with the `:left` command, so for the entire file, it'd be `:%left`
* `=i{` - indents based on curly brackets (syntax indenting will do this automatically)
* `<<` or `>>` - shift the current line backwards or forwards one indentation
* `<` or `>` - shift the current line backwards of forwards one indentation (Visual mode)

## Folding
* `:set foldmethod=indent` - folds based on indenting (I have this in my .vimrc)
* `zc` - closes fold cursor is in
* `zo` - opens fold cursor is in
* `zM` - closes all folds
* `zR` - opens all folds

## Sessions
* `:mks` - make a session
* `:mks!` - overwrite a session
* use vim-obsession + vim-prosession for automatic saving of sessions

## Tabs (USE BUFFERS, NOT TABS)
* `<C-w>T` - make focused split into a new tab
* `{i}gt` - move to the ith tab
* `:tabm {+/-}{i}` - move current tab +/- **i** tab spaces
* `gg=G` - indents an entire file based on your vimrc tab commands

## Searching
* `\#` - searches for the variable currently under the cursor
* `f{i}` - find and move cursor to the next character, **i**, on the given line. `F` finds backwards.
* `t{i}` - find and move the cursor to the character before the given one **i**. `T` does the same backwards.

## Various Other Commands
### Command Mode
* `:wa` - saves all tabs
* `:wqa` - saves all tabs and quits
* `:bd {x}` - deletes buffer **x**, where **x** is the file name or buffer \#

### Normal Mode
* `^` - go to beginning of line
* `$` - go to end of line
* `%` - when cursor is on opening bracket, goes to corresponding closing bracket
* `df{i}` - replace x with any symbol to delete until that symbol is found
* `D` - deletes to the end of the line (same as `d$`)
* `I` - enter insert mode at beginning of current line, disregarding leading whitespace
* `A` - enter insert mode at the end of current line

[This](http://vimregex.com/#pattern) is a good description of regex commands in vim

---

# tmux commands
* All commands prefixed by "tmux" are run from outside tmux, unless otherwise preceded by a colon
  * `<leader> - <C-b>` in fresh tmux installations. I use `<C-a>` and unbind `<C-b>` due to interactions with vim
  * `<leader>:` - run typed commands from within a tmux session

## Sessions
* `tmux new -s {session-name} `- start a new session from outside of tmux with **session-name** as the name
* `tmux kill-session -t {session-name}` - kills the session with **session-name** as the name
* `:tmux kill-session -a` - kills all session except current session (run from inside tmux session)
* `<leader> $` - rename a session
* `<leader> d` - detach and return to terminal
* `<leader> s` - list all sessions (from within tmux). same as `tmux ls` from inside terminal
* `tmux a` - attach last used session
* `tmux a -t {session-name}` - attach session with **session-name** as the name
* `<leader> (` or `<leader> )` - move to next or previous session, respectively

## Windows
* `<leader> c` - create new window
* `<leader> ,` - rename current window
* `<leader> &` - close current window
* `<leader> n` or `<leader> p` - move to next or previous window, respectively. I use `<leader> <C-h>` and `<leader> <C-l>`
* `:swap-window -s {i} -t {j}` - swap the **i**th window with the **j**th window
* `:swap-window -s {i} -{x}` - move the **i**th window **x** places

## Panes
* `<leader> %` or `<leader> "` - split window into panes vertically or horizontally. I use `<leader> |` or `<leader> -`
* `<leader> {` or `<leader> }` - move current pane left or right
* `<leader> z` - maximize or minimize current pane (toggle)
* `<leader> !` - convert pane into a new window
* `:joinp -s {i}.{j} -t {x}.{y}` - join the **j**th pane of window **i** with the **y**th pane of window **x**. Must always have a pane number, even if the window only has one pane (i.e. `:joinp -s 1.1 -t 2.1`)
* `<leader> ↑ or ↓ or ← or →` - resize window up, down, left, or right, respectively. can hold down <Ctrl> to iteratively resize
* `<leader> x` - close current pane
* `<leader> C-x` - toggle send command to all panes in window (this is set in my .tmux.conf file, bind C-x setw synchronize-panes)

## Buffer
* `:setw -g mode-keys vi` - use vi key-bindings in buffer (obviously this is in my .tmux.conf file)
* `<leader> [` - enter buffer mode
* `<leader> g` - goto top
* `<leader> G` - goto bottom
* `<leader> q` or `<leader> <Esc>` - quit buffer 

---

# Shell Commands
## Generic
### FZF
* `fzf | pbcopy` - fuzzy search with fzf, and chosen path is copied with `<ENTER>`
* `find (path) | fzf` - fuzzy search with fzf in a (path) directory you're not currently in. can optionally add `| pbcopy`
* `fd . (path) | fzf` - same as above, but with fd instead of find. can optionally add `| pbcopy`

## fish shell commands
### vim-key-bindings
#### vim commands which do NOT work
* `df<x>`
* search

### brace expansion
* stuff in curly braces will be expanded so each element becomes a new parameter
  * `echo input.{c,h,txt}` - outputs 'input.c input.h input.txt'
  * `mv *.{dat,txt} dir` - moves all files with .dat and .txt file types to directory
* I have a function to batch change filetypes of files (fish_scripts/rename.fish)

```fish
function rename
    set __bash_args $argv
        for file in *.$__bash_args[1]
        mv -v -- "$file" (basename $file .$__bash_args[1]).$__bash_args[2]
    end
end
```
