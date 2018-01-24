sed -i 's/define __[A-Z]*__/define __FULL__/g' reaxc_forces_sunway.cpp
SKIPDEP=1 make -j 32 sunway && bsub -debug -I -b -m 1 -p -share_size 6144 -host_stack 256 -priv_size 4 -cgsp 64 -n 1 ./lmp_sunway -in SUNWAY-TESTS/in.sunway.reaxc -var cwd SUNWAY-TESTS/ -sf sunway -var S 1 -var T 1 > interaction_list_full
sed -i 's/define __[A-Z]*__/define __HALF__/g' reaxc_forces_sunway.cpp
SKIPDEP=1 make -j 32 sunway && bsub -debug -I -b -m 1 -p -share_size 6144 -host_stack 256 -priv_size 4 -cgsp 64 -n 1 ./lmp_sunway -in SUNWAY-TESTS/in.sunway.reaxc -var cwd SUNWAY-TESTS/ -sf sunway -var S 1 -var T 1 > interaction_list_half

sed -i 's/\r//g' interaction_list_*
vimdiff interaction_list_*



