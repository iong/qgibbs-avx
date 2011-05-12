#!/bin/bash

print_cfg () {
cat << EOF
250 250 0.24 0.45 1.3 20
250 250 0.123 0.505 1.25 10
250 250 0.076 0.606 1.15 50
150 350 0.0281 0.698 1.00 150
150 350 0.0115 0.753 0.9 400
50 450 0.0028 0.820 0.75 2000
EOF
}

save_data() {
	d=$1
	rm -rf $d
	mkdir $d
	cat Um_*.dat > $d/Um.dat
	mv *.log $d
}
#echo "#!/bin/bash" > tmp_script
#print_cfg | awk '{print "./gibbs3h " $1 " " $2 " " $3 " " $4 " " $5 " 0 >" $5 ".log &"}' >> tmp_script
#echo "wait" >> tmp_script
#chmod 755 tmp_script
#./tmp_script
#save_data disp+vol@dest

echo "#!/bin/bash" > tmp_script
print_cfg | awk '{print "./gibbs3 " $1 " " $2 " " $3 " " $4 " " $5 " " $6 " >" $5 ".log &"}' >> tmp_script
echo "wait" >> tmp_script
chmod 755 tmp_script
./tmp_script
save_data disp+vol+swap@dest
