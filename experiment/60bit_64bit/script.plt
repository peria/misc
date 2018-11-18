set terminal png
set logscale x
set logscale y
set key left top

set output 'error.png'
set ylabel "Rounding error"
plot 'data.txt' index 0 using 1:6 ti "64bit 2+", \
     'data.txt' index 1 using 1:6 ti "60bit 2+", \
     'data.txt' index 2 using 1:6 ti "60bit 3+", \
     'data.txt' index 3 using 1:6 ti "60bit 5+", \
     'data.txt' index 4 using 1:6 ti "64bit 2-", \
     'data.txt' index 5 using 1:6 ti "60bit 2-", \
     'data.txt' index 6 using 1:6 ti "60bit 3-", \
     'data.txt' index 7 using 1:6 ti "60bit 5-"

set output 'error2.png'
set ylabel "Rounding error"
plot 'data.txt' index 0 using 1:6 ti "64bit +", \
     'data.txt' index 1:3 using 1:6 ti "60bit +", \
     'data.txt' index 4 using 1:6 ti "64bit -", \
     'data.txt' index 5:7 using 1:6 ti "60bit -"

set output 'time.png'
set ylabel "Time [ms]"
set xrange [1e+4:1e+9]
plot 'data.txt' index 0 using 1:2 ti "64bit 2+", \
     'data.txt' index 1 using 1:2 ti "60bit 2+", \
     'data.txt' index 2 using 1:2 ti "60bit 3+", \
     'data.txt' index 3 using 1:2 ti "60bit 5+", \
     'data.txt' index 4 using 1:2 ti "64bit 2-", \
     'data.txt' index 5 using 1:2 ti "60bit 2-", \
     'data.txt' index 6 using 1:2 ti "60bit 3-", \
     'data.txt' index 7 using 1:2 ti "60bit 5-"
