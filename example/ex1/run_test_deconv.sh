
rm -f wave_at_depth_*.txt
rm -f test_deconvolution.out

make test_deconvolution
./test_deconvolution.out

python plot_deconvolution.py wave_at_depth_*_acc.txt


