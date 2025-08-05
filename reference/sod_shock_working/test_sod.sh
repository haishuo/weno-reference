#!/bin/bash
echo "Testing Sod Shock Tube Implementation..."

# Compile
gfortran -O3 -march=native -fdefault-real-8 -fdefault-double-8 \
         -fno-trapping-math -ffpe-summary=none -o weno weno_with_sod.f

# Run Sod test
printf "3\n20 10\n0.3\n10\n0.1\n4\n0\n" | ./weno

# Verify results
echo "Verification:"
echo "Left state (x=4.5):"
awk 'NR>1 && $1 == 4.5 {print "  x=" $1 ", rho=" $3 ", p=" $4}' fort.9 | head -1

echo "Right state (x=5.0):"
awk 'NR>1 && $1 == 5.0 {print "  x=" $1 ", rho=" $3 ", p=" $4}' fort.9 | head -1

echo "Expected: Left rho=1.0 p=1.0, Right rho=0.125 p=0.1"
