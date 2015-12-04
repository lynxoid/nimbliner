EXE=mapper
rm $EXE
make $EXE
/usr/bin/time -v ./$EXE query 20 $1 $2 $3
