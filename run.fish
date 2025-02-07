#! /usr/bin/env fish
set default "(x^2+y^2)^5-(x^4-6x^2y^2+y^4)^2
xy(x^2-y^2)"

make clean && make -j8 && 
if test "$argv[1]" = 'dv' -o "$argv[1]" = 'vd'
    echo $default | valgrind ./prog
else if test "$argv[1]" = 'v'
    valgrind ./prog
else if test "$argv[1]" = 'd'
    echo $default | ./prog
else
    ./prog
end &&
ffmpeg -hide_banner -loglevel warning -framerate 20 -pattern_type glob -i '*.bmp' -r 15 -vf scale=1440:-1 -c:v libx264 -pix_fmt yuv420p -y out.mp4
rm *.bmp
