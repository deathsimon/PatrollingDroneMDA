# Check if input file is provided
if (strlen(input_file) == 0) {
    print "Usage: gnuplot -e \"input_file='data.txt'\" plot_script.gp"
    exit
}

# Set the terminal to png
set terminal pngcairo dashed size 800,600 enhanced font 'Verdana,16' 

# Set the output file
set output 'result.png'

# Set the title of the graph
# set title "Comparison of Greedy, TSP, and Enforced Algorithms"

# Set the labels for the axes
set xlabel "Scenario ID"
set ylabel "Normalized Maximum Age-of-Information"
set xtics 20

set yrange [0.98:1.5]

# Set the key (legend)
set key top right

# Define custom styles for lines
set style line 1 lc rgb '#9400d3' lt 1 lw 3 # Red line, solid, width 3
set style line 2 lc rgb '#e69f00' lt 1 lw 3 # Green line, solid, width 3
set style line 3 lc rgb '#009e73' lt 1 lw 3 # Blue line, solid, width 3

# Plot the data
plot input_file using 1 with lines linestyle 3 title 'Greedy', \
     '' using 2 with lines linestyle 2 title 'SRTT', \
     '' using 3 with lines linestyle 1 title 'Enforced'