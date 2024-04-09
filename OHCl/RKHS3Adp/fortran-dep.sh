#!/bin/bash

for src in $*; do
    ext=${src##*.}
    base=${src%.*}
    if [ "$base.$ext" == "$src" ]; then
	echo "\$(O)/$base.o: $src"
    fi
done
