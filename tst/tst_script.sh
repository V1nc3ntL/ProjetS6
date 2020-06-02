#!/bin/sh
sz=$1
res_f="matlab_res.txt"
sig="signal.txt"
# Passer en argument à matlab et récupérer les valeurs
matlab_test()
{
chmod +770 $1;
smp=$(tr '\n[0-9]' ',[0-9]' < "$1") 
#On enlève le dernier caractère
smp=${smp%?}
#echo $smp
str="n=fft([$smp]);fprintf(\"%f\n%f\n\",n,imag(n))"
matlab -batch "$str" > $res_f
}

./TEST init "$sz" ;
matlab_test $sig 
./TEST $sig $res_f ;
