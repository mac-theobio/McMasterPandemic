#!/bin/csh -f

if ($#argv == 0) then
   set ff = `ls`
else
   set ff = "$argv"
endif

echo files to be changed:  $ff

foreach f ($ff)
    echo $f
    mv $f $f.1
    sed 's/gamma/sigma/g' $f.1 > $f.2
    sed 's/lambda/gamma/g' $f.2 > $f.3
    mv -f $f.3 $f
    rm $f.2
    git diff -U0 --word-diff --no-index $f.1 $f
    echo ------------------------------------------------
end
