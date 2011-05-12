#!/bin/bash

print_cfg () {
cat << EOF
250 250 0.35 0.35 1.3 20
250 250 0.3 0.3 1.25 10
250 250 0.2 0.2 1.15 50
150 350 0.06 0.60 1.00 150
150 350 0.06 0.60 0.9 400
50 450 0.001 0.6 0.75 2000
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
exit
./tmp_script
save_data disp+vol+swap@dest
