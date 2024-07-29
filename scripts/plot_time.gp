# Check if input file is provided
if (strlen(input_file) == 0) {
    print "Usage: gnuplot -e \"input_file='data.txt'\" plot_script.gp"
    exit
}

# Set the terminal to png
set terminal pngcairo size 1000,600 enhanced font 'Verdana,16'

# Set the output file
set output 'bar_chart.png'

# Set the title of the graph
# set title "Comparison of Records"

# Set the labels for the axes
# set xlabel "Records"
set ylabel "Execution Time (ms)"

# Set log scale for y-axis
set logscale y

# Set the style for the bars
set style data histograms
set style histogram clustered gap 1
set style fill solid border -1

# Set the xtics
#set xtics rotate by -30

set key top left

# Read the data file
plot input_file using 2:xtic(1) title 'Greedy' linecolor rgb "#009e73" fs pattern 3, \
     '' using 3 title 'SRTT' linecolor rgb "#e69f00" fs pattern 4, \
     '' using 4 title 'Enforced' linecolor rgb "#9400d3" fs pattern 1, \
     '' using 5 title 'DP' linecolor rgb "#56b4e9" fs pattern 2
