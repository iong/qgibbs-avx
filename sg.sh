#!/bin/bash

print_cfg () {
cat << EOF
250 250 1.3e-2 1.3e-2 44 100
250 250 1.2e-2 1.3e-2 42 100
250 250 1e-2 1e-2 40 100
250 250 7e-3 7e-3 37 500
150 350 2e-3 0.025 32 1500
150 350 2e-3 0.025 29 4000
50 450 3.5e-5 0.025 24 4000
EOF
}

save_data() {
	d=$1
	rm -rf $d
	mkdir $d
	cat Um_*.dat > $d/Um.dat
	mv Um_*.dat *.log $d
}

echo "#!/bin/bash" > tmp_script
print_cfg | sed -e 's/^/.\/gibbs3 /' -e 's/$/ \&/' >> tmp_script
echo "wait" >> tmp_script
chmod 755 tmp_script
./tmp_script
save_data SG
