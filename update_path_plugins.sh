#!/usr/local/bin/fish

set cwd (pwd)
cd ~/.vim/bundle/

for i in */
    cd "$i"
    git pull
    cd -
end

cd $cwd
set -e cwd
