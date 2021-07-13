#!/bin/bash
# location
echo $1
cd $1
echo $(pwd)

git pull
git config --global user.name 'calibration-bot'
git config --global user.email 'calibration-bot@users.noreply.github.com'
git add --all
git commit -am "Weekly regression testing calibration [skip ci]"
git push