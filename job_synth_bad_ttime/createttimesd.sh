#! /bin/bash
#
# Create traveltimes using the taup package

echo "=> Quake coordinates:"
echo "Depth: 5.0 km"
echo "Quake: 10,10 [km]"
echo "Quake: 0.0899,0.0899 [deg]"
echo
echo "=> Station coordinates:"
echo "GAK2:   0, 0 [deg]"
echo "GAK3:   0, 0.045 [deg]"
echo "GAK4:   0.045, 0.045 [deg]"
echo
echo "=>  GAK2"

depth=5

taup_time -mod regional -h ${depth} -deg 0.1271 -ph p,s3p,p^0Pv3p,p^0Pv3p^0Pv3p

echo "=>  GAK3"
taup_time -mod regional -h ${depth} -deg 0.1005 -ph p,s3p,p^0Pv3p,p^0Pv3p^0Pv3p

echo "=>  GAK4"
taup_time -mod regional -h ${depth} -deg 0.0635 -ph p,s3p,p^0Pv3p,p^0Pv3p^0Pv3p




