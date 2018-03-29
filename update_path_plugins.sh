#!/usr/local/bin/fish
# short script to update all vim plugins (tested for pathogen)
# written for fish shell

set cwd (pwd)
cd ~/.vim/bundle/

if test -e ~/Google\ Drive/github_repos/mjv_personal_code/installed_vim_plugins.txt
    rm ~/Google\ Drive/github_repos/mjv_personal_code/installed_vim_plugins.txt
end

for i in */
    echo "$i" >> ~/Google\ Drive/github_repos/mjv_personal_code/installed_vim_plugins.txt
    cd "$i"
    git pull
    cd -
end

echo "installed_vim_plugins.txt written to personal GitHub repo"
cd $cwd
set -e cwd
