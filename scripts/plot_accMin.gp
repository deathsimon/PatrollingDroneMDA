# Check if input file is provided
if (strlen(input_file) == 0) {
    print "Usage: gnuplot -e \"input_file='data.txt'\" plot_script.gp"
    exit
}

# Set the terminal to PNG
set terminal pngcairo size 1000,600 enhanced font 'Verdana,16'

# Set the output file
set output 'stacked_bar_chart.png'

# Set the title of the graph
# set title "Stacked Bar Chart"

# Set the labels for the axes
# set xlabel "Records"
set ylabel "\# of Times"

# Set the grid
set grid ytics

# Set the style for the bars
set style data histograms
set style histogram rowstacked gap 100
set boxwidth 0.75 relative
set style fill solid border -1

# Set the xtics to use the first column (record names) as labels
# set xtics rotate by -45

# Set the legend
set key above horizontal

# Read the data file
plot input_file using 8:xtic(1) title 'Enforced' linecolor rgb "#9400d3" fs pattern 1, \
     '' using 7 title 'SRTT' linecolor rgb "#e69f00" fs pattern 4, \
     '' using 6 title 'Greedy' linecolor rgb "#009e73" fs pattern 3
