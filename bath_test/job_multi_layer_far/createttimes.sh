#! /bin/bash
#
# Create traveltimes using the taup package

seaf=4.3 # depth of seafloor interface

echo "=> Quake coordinates:"
echo "Depth: 5.0 km"
echo "Quake: 10,10 [km]"
echo
echo "=> Station coordinates:"
echo "GAK2:   0, 0 [km]"
echo "GAK3:   0, 5 [km]"
echo "GAK4:   5, 5 [km]"
echo
echo "=>  GAK2"
taup_time -mod regional -h 5.0 -km 14.142 -ph "p,s${seaf}p,p^0Pv${seaf}p,p^0Pv${seaf}p^0Pv${seaf}p"

echo "=>  GAK3"
taup_time -mod regional -h 5.0 -km 11.180 -ph "p,s${seaf}p,p^0Pv${seaf}p,p^0Pv${seaf}p^0Pv${seaf}p"

echo "=>  GAK4"
taup_time -mod regional -h 5.0 -km 7.071 -ph "p,s${seaf}p,p^0Pv${seaf}p,p^0Pv${seaf}p^0Pv${seaf}p"




