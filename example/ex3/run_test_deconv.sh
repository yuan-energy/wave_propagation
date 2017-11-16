
rm -f wave_at_depth_*.txt
rm -f test_deconvolution.out

make test_deconvolution
./test_deconvolution.out

python plot_deconvolution.py 
python plot.py  wave_at_depth_0_acc.txt
python plot.py  wave_at_depth_250_acc.txt


