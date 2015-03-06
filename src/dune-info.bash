#!/bin/bash

exec &> dune-info.dat

pushd . >/dev/null

# Find base of git directory
while [ ! -d .git ] && [ ! `pwd` = "/" ]; do cd ..; done

# Go to Dune root dir
cd ..

# Collect information about all Dune module using dunecontrol
#echo "== SVN info: "
#dune-common/bin/dunecontrol svn info
#echo

#echo "== SVN diff: "
#dune-common/bin/dunecontrol svn diff
#echo

echo "== Git info: "
dune_dirs=`find . -maxdepth 1 -type d -name "dune-*"`
echo $dune_dirs
for dir in $dune_dirs ; do
 echo "--- calling git log for $dir ---"
 cd $dir
 git log --max-count=1
 cd ..
 echo "--- $dir done ---"
done
echo

echo "== Git diff: "
for dir in $dune_dirs ; do
 echo "--- calling git diff for $dir ---"
 cd $dir
 git diff HEAD
 cd ..
 echo "--- $dir done ---"
done
echo

popd >/dev/null
