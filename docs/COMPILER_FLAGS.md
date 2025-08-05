# Critical FORTRAN Compiler Flags

## ⚠️ ESSENTIAL FLAGS

The WENO FORTRAN code **REQUIRES** these specific compiler flags to work correctly:

```bash
gfortran -O3 -march=native \
         -fdefault-real-8 -fdefault-double-8 \
         -fno-trapping-math -ffpe-summary=none \
         -o weno weno.f
```

## Why Each Flag is Needed

| Flag | Purpose | Critical? |
|------|---------|-----------|
| `-fdefault-real-8` | Forces all REAL variables to double precision | **ESSENTIAL** |
| `-fdefault-double-8` | Forces all DOUBLE PRECISION variables to 8 bytes | **ESSENTIAL** |
| `-fno-trapping-math` | Allows IEEE floating-point operations without trapping | **REQUIRED** |
| `-ffpe-summary=none` | Suppresses floating-point exception summaries | **REQUIRED** |
| `-O3` | Optimization level 3 for performance | Recommended |
| `-march=native` | Optimize for current CPU architecture | Recommended |

## What Happens Without These Flags

**Missing `-fdefault-real-8 -fdefault-double-8`:**
- ❌ **NaN values everywhere** in output
- ❌ **Numerical instability** in WENO reconstruction
- ❌ **Complete failure** of shock-capturing

**Missing `-fno-trapping-math -ffpe-summary=none`:**
- ❌ **Floating-point exceptions** and crashes
- ❌ **IEEE_INVALID_FLAG** errors
- ❌ **Premature termination** with SIGFPE

## Historical Context

Chi-Wang Shu's original WENO code (1990s) was written for older FORTRAN compilers that:
- Had different default precision behavior
- Were more permissive with floating-point operations
- Used different optimization strategies

Modern gfortran requires explicit flags to maintain compatibility.

## Test Commands

### Compilation Test
```bash
cd fortran/modified
gfortran -O3 -march=native -fdefault-real-8 -fdefault-double-8 \
         -fno-trapping-math -ffpe-summary=none -o weno weno.f
```

### Validation Test (Sod Shock Tube)
```bash
printf "3\n20 10\n0.3\n10\n0.1\n4\n0\n" | ./weno
```

**Expected Output:**
- Clean execution without errors
- `fort.9` shows density jump from 1.0 to 0.125 at x=5.0
- No NaN values in any output files

### Quick Verification
```bash
# Check for successful shock discontinuity
awk 'NR>1 && $1 < 5.0 {left=$3} NR>1 && $1 >= 5.0 && !right {right=$3} 
     END {printf "Left: %.3f, Right: %.3f, Ratio: %.1f\n", left, right, left/right}' fort.9
```

Should output: `Left: 1.000, Right: 0.125, Ratio: 8.0`

## Alternative Compiler Approaches

### For Debugging
```bash
# Add debug information and stricter checking
gfortran -g -fcheck=all -fdefault-real-8 -fdefault-double-8 \
         -fno-trapping-math -ffpe-summary=none -o weno_debug weno.f
```

### For Maximum Performance
```bash
# Additional optimization flags
gfortran -O3 -march=native -mtune=native -funroll-loops \
         -fdefault-real-8 -fdefault-double-8 \
         -fno-trapping-math -ffpe-summary=none -o weno_fast weno.f
```

### For Legacy Systems
```bash
# If march=native causes issues
gfortran -O3 -fdefault-real-8 -fdefault-double-8 \
         -fno-trapping-math -ffpe-summary=none -o weno weno.f
```

## Troubleshooting

### Problem: "gfortran: command not found"
**Solution:** Install gfortran
```bash
# Ubuntu/Debian
sudo apt install gfortran

# macOS with Homebrew  
brew install gcc

# Red Hat/CentOS
sudo yum install gcc-gfortran
```

### Problem: Still getting NaN despite correct flags
**Solution:** Check that you're using the modified version
```bash
# Verify you have the shock tube modifications
grep -n "Override for shock problems" weno.f
```

### Problem: Warnings about "Fortran 2018 deleted feature"
**Solution:** These are harmless warnings about old-style DO loops
```
Warning: Fortran 2018 deleted feature: Shared DO termination label
```
These warnings can be ignored - they don't affect functionality.

## Summary

**Minimum Required Flags:**
```bash
gfortran -fdefault-real-8 -fdefault-double-8 -fno-trapping-math -ffpe-summary=none -o weno weno.f
```

**Recommended Full Command:**
```bash
gfortran -O3 -march=native -fdefault-real-8 -fdefault-double-8 -fno-trapping-math -ffpe-summary=none -o weno weno.f
```

**Remember:** These flags are the result of extensive debugging and are **absolutely critical** for the code to work correctly!
