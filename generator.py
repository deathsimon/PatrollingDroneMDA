import json
import os
import random
import math

def generate_coordinates_grid(nodes, sideLength):
    coordinates = []
    step = sideLength / 4

    for i in range(4):
        for j in range(4):
            for _ in range(nodes // 16 +1):            
                x = math.ceil(random.uniform((-1)*sideLength/2 + i * step, (-1)*sideLength/2 + (i + 1) * step))
                y = math.ceil(random.uniform((-1)*sideLength/2 + j * step, (-1)*sideLength/2 + (j + 1) * step))
                coordinates.append({"x": x, "y": y})
    return random.sample(coordinates,nodes)

def generate_coordinates_cluster(nodes, sideLength):
    coordinates = []
    step = sideLength / 4

    if nodes <= 8:
        clusters = 1
    else:
        clusters = 4

    chosen_areas = random.sample(range(16), clusters)
    for area in chosen_areas:
        for _ in range(nodes // clusters +1):
            i = area // 4
            j = area % 4
            x = math.ceil(random.uniform((-1)*sideLength/2 + i * step, (-1)*sideLength/2 + (i + 1) * step))
            y = math.ceil(random.uniform((-1)*sideLength/2 + j * step, (-1)*sideLength/2 + (j + 1) * step))
            coordinates.append({"x": x, "y": y})

    # for k in range(nodes % 4):
    #     i = chosen_areas[k] // 4
    #     j = chosen_areas[k] % 4
    #     x = random.randrange((-1)*sideLength/2 + i * step, (-1)*sideLength/2 + (i + 1) * step)
    #     y = random.randrange((-1)*sideLength/2 + j * step, (-1)*sideLength/2 + (j + 1) * step)            
    #     coordinates.append({"x": x, "y": y})
    return random.sample(coordinates,nodes)

def generate_coordinates_outlier(nodes, sideLength):
    coordinates = []
    step = sideLength / 2
    
    chosen_areas = random.sample(range(4), 2)
    outlier_area = chosen_areas[0]
    cluster_area = chosen_areas[1]

    i_outlier = outlier_area // 2
    j_outlier = outlier_area % 2
    x_outlier = random.randrange((-1)*sideLength/2 + i_outlier * step, (-1)*sideLength/2 + (i_outlier + 1) * step)
    y_outlier = random.randrange((-1)*sideLength/2 + j_outlier * step, (-1)*sideLength/2 + (j_outlier + 1) * step)
    coordinates.append({"x": x_outlier, "y": y_outlier})

    i_cluster = cluster_area // 2
    j_cluster = cluster_area % 2
    for _ in range(nodes - 1):
        x = random.randrange((-1)*sideLength/2 + i_cluster * step, (-1)*sideLength/2 + (i_cluster + 1) * step)
        y = random.randrange((-1)*sideLength/2 + j_cluster * step, (-1)*sideLength/2 + (j_cluster + 1) * step)
        coordinates.append({"x": x, "y": y})
    return coordinates

def main():
    # Read input JSON file
    with open("config.json", "r") as input_file:
        input_data = json.load(input_file)

    mode = input_data["mode"]
    nodes = input_data["data nodes"]
    cases = input_data["cases"]
    sideLength = input_data["length"]    

    # Generate output filename
    output_filename = f"{mode}-{nodes}-{cases}.json"

    # Check if output file already exists
    if os.path.exists(output_filename):
        overwrite = input(f"The file {output_filename} already exists. Do you want to overwrite it? (yes/no): ")
        if overwrite.lower() != "yes":
            print("Operation aborted.")
            return

    output_data = {"cases": cases, "scenarios": []}

    # Generate coordinates
    for _ in range(cases):
        coordinates = [{"x": 0, "y": 0}]
        if mode == "grid":
            coordinates.extend(generate_coordinates_grid(nodes, sideLength))
        elif mode == "cluster":
            coordinates.extend(generate_coordinates_cluster(nodes, sideLength))
        elif mode == "outlier":
            coordinates.extend(generate_coordinates_outlier(nodes, sideLength))
        else:
            print(f"Unknown mode: {mode}")
            return
        output_data["scenarios"].append({"N": nodes+1,"points": coordinates})

    # Write output JSON file
    with open(output_filename, "w") as output_file:
        json.dump(output_data, output_file, indent=4)

    print(f"Coordinates generated and saved to {output_filename}")

if __name__ == "__main__":
    main()