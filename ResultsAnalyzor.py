import argparse
import json
import numpy as np
#import matplotlib.pyplot as plt

# Function to read JSON file
def read_json_file(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
    return data

# Function to normalize data
def normalize_data(data):
    normalized_data = []
    for record in data:
        optimal = record['Optimal']['MDA']
        normalized_record = {
            'Greedy': record['Greedy']['MDA'] / optimal,
            'MST': record['MST']['MDA'] / optimal,
            'Enforced': record['Enforced']['MDA'] / optimal,
            'Optimal': 1  # Normalized to itself
        }
        normalized_data.append(normalized_record)
    return normalized_data

# Function to sort data based on normalized Enforced value
def sort_data_by_enforced(normalized_data):
    sorted_data = sorted(normalized_data, key=lambda x: x['Enforced'])
    return sorted_data

# Function to calculate mean and variance
def calculate_statistics(normalized_data):
    greedy_values = [record['Greedy'] for record in normalized_data]
    mst_values = [record['MST'] for record in normalized_data]
    enforced_values = [record['Enforced'] for record in normalized_data]
    
    statistics = {
        'Greedy': {'mean': np.mean(greedy_values), 'variance': np.var(greedy_values)},
        'MST': {'mean': np.mean(mst_values), 'variance': np.var(mst_values)},
        'Enforced': {'mean': np.mean(enforced_values), 'variance': np.var(enforced_values)}
    }
    
    return statistics

# Function to write normalized data to a file
def write_normalized_data_to_file(normalized_data, output_file):
    with open(output_file, 'w') as file:
        line = f"{'Greedy'}\t{'MST'}\t{'Enforced'}\n"
        file.write(line)
        for record in normalized_data:
            line = f"{record['Greedy']}\t{record['MST']}\t{record['Enforced']}\n"
            file.write(line)

def avg_time(data):
    avg_time = []

    sum_time_greedy = 0.0
    sum_time_mst = 0.0
    sum_time_enforced = 0.0
    sum_time_optimal = 0.0
    
    for record in data:
        sum_time_greedy += record['Greedy']['Time']
        sum_time_mst += record['MST']['Time']
        sum_time_enforced += record['Enforced']['Time']
        sum_time_optimal += record['Optimal']['Time']

    avg_time.append(sum_time_greedy/len(data))
    avg_time.append(sum_time_mst/len(data))
    avg_time.append(sum_time_enforced/len(data))
    avg_time.append(sum_time_optimal/len(data))
        
    return avg_time

def count_best(data):
    best={0,0,0}
    for record in data:
        if record['Greedy']['MDA'] < record['MST']['MDA'] and record['Greedy']['MDA'] < record['Enforced']['MDA']:
            best[0] += 1
        elif record['MST']['MDA'] < record['Greedy']['MDA'] and record['MST']['MDA'] < record['Enforced']['MDA']:
            best[1] += 1
        else:
            best[2] += 1

    return best
    

# Function to plot the data
# def plot_data(normalized_data):
#     indices = range(len(normalized_data))
#     greedy_values = [record['Greedy'] for record in normalized_data]
#     mst_values = [record['MST'] for record in normalized_data]
#     enforced_values = [record['Enforced'] for record in normalized_data]

#     plt.figure(figsize=(10, 6))
#     plt.plot(indices, greedy_values, label='Greedy', marker='o')
#     plt.plot(indices, mst_values, label='MST', marker='o')
#     plt.plot(indices, enforced_values, label='Enforced', marker='o')
    
#     plt.xlabel('Record Index')
#     plt.ylabel('Normalized Value')
#     plt.title('Normalized Values of Greedy, MST, and Enforced to Optimal')
#     plt.legend()
#     plt.grid(True)
#     plt.show()

# Main function
def main():
    parser = argparse.ArgumentParser(description='Process and output normalized data from a JSON file.')
    parser.add_argument('input_file', type=str, help='Path to the input JSON file')
    parser.add_argument('output_file', type=str, help='Path to the output file to save normalized data')
    
    args = parser.parse_args()
    input_file = args.input_file
    output_file = args.output_file

    #file_path = 'output-outlier-8-100.json'
    data = read_json_file(input_file)
    normalized_data = normalize_data(data)
    sorted_data = sort_data_by_enforced(normalized_data)
    statistics = calculate_statistics(sorted_data)
    
    print("Mean and Variance of Normalized Values:")
    for key, stats in statistics.items():
        print(f"{key} - Mean: {stats['mean']}, Variance: {stats['variance']}")
    
    #plot_data(sorted_data)
    write_normalized_data_to_file(sorted_data, output_file)

    timeResults = avg_time(data)
    print("Average Time:")
    print(f"Greedy: {timeResults[0]}")
    print(f"MST: {timeResults[1]}")
    print(f"Enforced: {timeResults[2]}")
    print(f"Optimal: {timeResults[3]}")

    best = count_best(data)
    print("Best Results:")
    print(f"Greedy: {best[0]}")
    print(f"MST: {best[1]}")
    print(f"Enforced: {best[2]}")

if __name__ == "__main__":
    main()