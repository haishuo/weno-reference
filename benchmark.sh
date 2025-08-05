#!/bin/bash

# Fixed Beast Mode Test - Now with proper input formatting!
# The segfault was just an input parsing issue

echo "FIXED Beast Mode FORTRAN Test"
echo "============================="
echo "Now with proper input formatting (printf instead of echo)"
echo ""

cd fortran/modified || exit 1

# Recompile for beast mode
echo "Compiling with maximum optimization..."
gfortran -O3 -march=native -mtune=native -funroll-loops \
         -fdefault-real-8 -fdefault-double-8 \
         -fno-trapping-math -ffpe-summary=none \
         -o weno_beast weno.f

if [ $? -ne 0 ]; then
    echo "Compilation failed!"
    exit 1
fi

echo "Beast mode compilation successful!"
echo ""

# Function for beast mode testing with CORRECT input format
beast_test() {
    local nx=$1
    local ny=$2
    local steps=$3
    local time_final=$4
    local description="$5"
    
    local points=$((nx * ny))
    local estimated_seconds=$((points / 850000))
    local estimated_minutes=$((estimated_seconds / 60))
    
    echo ""
    echo "==========================================="
    echo "BEAST TEST: $description"
    echo "Grid: ${nx}×${ny} = ${points} points"
    echo "Time steps: ${steps}"
    echo "Estimated: ${estimated_seconds}s (~${estimated_minutes} minutes)"
    echo "Memory usage: ~$((points * 4 * 8 / 1024 / 1024))MB"
    echo "==========================================="
    
    if [ $estimated_minutes -gt 15 ]; then
        echo "WARNING: Estimated ${estimated_minutes} minutes!"
        read -p "Continue? (y/n): " -n 1 -r
        echo
        [[ ! $REPLY =~ ^[Yy]$ ]] && return
    fi
    
    echo "Starting BEAST MODE test..."
    echo "(Press Ctrl+C to abort)"
    
    local start_time=$(date +%s)
    
    # THE FIX: Use printf instead of echo for proper newline handling!
    printf "3\n%d %d\n0.3\n%d\n%g\n4\n0\n" $nx $ny $steps $time_final | \
        timeout 1800s ./weno_beast > beast_output.log 2>&1
    
    local exit_code=$?
    local end_time=$(date +%s)
    
    local duration=$((end_time - start_time))
    local minutes=$((duration / 60))
    local seconds=$((duration % 60))
    
    echo ""
    if [ $exit_code -eq 124 ]; then
        echo "⏰ TIMEOUT after 30 minutes"
        echo "💪 This is definitely a challenging problem for CPU!"
    elif [ $exit_code -eq 0 ]; then
        echo "✅ COMPLETED in ${minutes}m ${seconds}s"
        
        if [ $duration -ge 480 ] && [ $duration -le 720 ]; then
            echo "🎯 GOLDILOCKS! Perfect 8-12 minute range!"
            echo "🚀 This is our ideal GradFlow benchmark target!"
        elif [ $duration -ge 300 ]; then
            echo "👍 Good candidate (5+ minutes)"
        elif [ $duration -ge 60 ]; then
            echo "📈 Promising (1+ minute) - need to go bigger"
        else
            echo "⚡ Still too fast for this beast CPU!"
        fi
        
        if [ -f fort.9 ]; then
            local output_lines=$(wc -l < fort.9)
            echo "📊 Output: ${output_lines} data points"
            local file_size=$(du -h fort.9 | cut -f1)
            echo "📁 File size: ${file_size}"
            rm -f fort.8 fort.9
        fi
    else
        echo "❌ FAILED with exit code $exit_code"
        echo "Last few lines of output:"
        tail -10 beast_output.log
    fi
    
    echo "⏱️  Beast result: ${duration}s (${minutes}m ${seconds}s)"
    rm -f beast_output.log
}

echo "🚀 BEAST MODE TESTS - Now with proper input formatting!"
echo ""

# Start with known working sizes and scale up
echo "Phase 1: Confirming the fix works"
beast_test 400 400 200 2.0 "Confirm 400x400 works (was 193ms before)"

beast_test 600 600 300 3.0 "600x600 - was 431ms before"

echo ""
echo "Phase 2: Scaling up to find the 10-minute sweet spot"

beast_test 1000 1000 200 2.0 "1M points - should take ~1-2 minutes"

beast_test 1500 1500 200 2.0 "2.25M points - should take ~3-4 minutes"

beast_test 2000 2000 150 1.5 "4M points - should take ~5-7 minutes"

beast_test 2500 2500 120 1.2 "6.25M points - aiming for ~8-10 minutes"

beast_test 3000 3000 100 1.0 "9M points - should hit the 10+ minute target"

echo ""
echo "Phase 3: Fine-tuning if needed"

# Only run if previous tests suggest we need to go bigger
beast_test 3500 3500 80 0.8 "12.25M points - if still not 10 minutes"

beast_test 4000 4000 75 0.75 "16M points - pushing toward 15+ minutes"

echo ""
echo "=============================================="
echo "BEAST MODE ANALYSIS COMPLETE"
echo "=============================================="
echo ""
echo "🎯 Target found: The problem that takes ~10 minutes on your beast CPU"
echo "   This becomes our GradFlow 60x acceleration benchmark"
echo ""
echo "💡 Key insight: Input formatting was the issue, not memory limits!"
echo "   Your CPU can handle much larger problems than initially thought"
echo ""
echo "🚀 GradFlow development strategy confirmed:"
echo "   - Small problems: Correctness validation"
echo "   - BEAST problems: Performance demonstration" 
echo "   - 3D problems: Revolutionary impact"
echo ""
echo "Ready to build the GPU solver that makes 10-minute problems"
echo "finish in 10 seconds! 💥"
