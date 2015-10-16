#!/bin/sh
cp ./scripts/* ./mcc
cp ./functions/* ./mcc
matlab -nojvm -nodisplay -nosplash -r "compilesimcode"
