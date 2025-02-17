import pandas as pd
import numpy as np
import joblib
import sys
import os

# Ensure correct usage
if len(sys.argv) != 5:
    print("Usage: python ml_script2.py <ml_scripts_dir> <input_snp_file> <bed_file> <output_predictions>")
    sys.exit(1)

# Get input arguments
ml_scripts_dir = sys.argv[1]
input_snp_file = sys.argv[2]
bed_file = sys.argv[3]
output_predictions = sys.argv[4]

# Ensure the model file exists in the ML directory
model_path = os.path.join(ml_scripts_dir, "best_model_top30.pkl")
if not os.path.exists(model_path):
    print(f"Error: Model file '{model_path}' not found.")
    sys.exit(1)

# Load model from ML directory
model = joblib.load(model_path)

# Load BED file for filtering SNP positions
hotspot_bed = pd.read_csv(bed_file, sep="\t", header=None, names=["Chromosome", "Start", "End"])
hotspot_bed["Start"] = hotspot_bed["Start"] + 1  # Convert to 1-based indexing

# CGR Encoding Function
def get_final_cgr_position(snp_vector, grid_size=8):
    x, y = grid_size // 2, grid_size // 2
    for snp in snp_vector:
        if snp == 0:
            x, y = x // 2, y // 2
        elif snp == 1:
            x, y = (x + grid_size) // 2, (y + grid_size) // 2
        elif snp == 2:
            x, y = (x + grid_size) // 2, y // 2
    distance = np.sqrt((x - (grid_size // 2))**2 + (y - (grid_size // 2))**2)
    return distance

# Process input SNP file
def create_features_from_snp_file(snp_file, bed_data):
    print("Processing SNP file:", snp_file)
    snp_data = pd.read_csv(snp_file, sep="\t", header=None)
    snp_data.columns = ["Chromosome", "Position"] + [f"Sample_{i}" for i in range(snp_data.shape[1] - 2)]

    # Prepare features for each sample
    num_samples = snp_data.shape[1] - 2
    features_list = [[] for _ in range(num_samples)]

    # Match SNP positions with given BED file
    for _, row in bed_data.iterrows():
        chrom, start, end = row["Chromosome"], row["Start"], row["End"]
        match = snp_data[(snp_data["Chromosome"] == chrom) & (snp_data["Position"].between(start, end))]

        for i in range(num_samples):
            if not match.empty:
                snp_vector = match.iloc[:, i + 2].values
                cgr_feature = get_final_cgr_position(snp_vector)
                features_list[i].append(cgr_feature)
            else:
                features_list[i].append(0)

    return [np.array(features).reshape(1, -1) for features in features_list]  # Ensure model input format

# Make prediction for all samples
def predict_from_snp_file(snp_file, bed_file, output_file):
    feature_sets = create_features_from_snp_file(snp_file, hotspot_bed)
    predictions = []

    with open(output_file, "w") as out_f:
        out_f.write("Sample\tPrediction\tConfidence (%)\n")
        for i, features in enumerate(feature_sets):
            prediction = model.predict(features)
            probabilities = model.predict_proba(features)
            confidence = probabilities.max() * 100
            predictions.append((prediction[0], confidence))
            out_f.write(f"Sample_{i + 1}\t{prediction[0]}\t{confidence:.2f}\n")
            print(f"Predicted Cultivar for Sample {i + 1}: {prediction[0]} (Confidence: {confidence:.2f}%)")

    return predictions

# Run prediction
predict_from_snp_file(input_snp_file, bed_file, output_predictions)
