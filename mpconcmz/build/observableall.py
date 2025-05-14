import pandas as pd
import numpy as np
import os
import re

def extract_N_from_filename(filename):
    """Extract N value from filename like results20...csv → 20"""
    match = re.search(r'results(\d+)', filename)
    if match:
        return int(match.group(1))
    return None

def process_csv_file(input_path, output_path, N):
    """Process a single CSV file with N-dependent scaling and derivative calculation."""
    try:
        df = pd.read_csv(input_path)
        
        # Divide C1, C2, C3 and their fluctuations by N
        for col in ['C1', 'C2', 'C3']:
            if col in df.columns:
                df[col] = df[col] / N
            if f"{col}fluc" in df.columns:
                df[f"{col}fluc"] = df[f"{col}fluc"] / N
        
        # Calculate dC1/dT (assuming C1 represents <n_c>)
        if 'C1' in df.columns and 'Temperature' in df.columns:
            T = df['Temperature'].values
            C1 = df['C1'].values
            
            # Compute derivative using central differences
            dC1_dT = np.gradient(C1, T)
            df['dC1_dT'] = dC1_dT
            
            # Optional: Smooth the derivative if data is noisy
            # from scipy.interpolate import UnivariateSpline
            # spline = UnivariateSpline(T, C1, s=0.5)
            # df['dC1_dT_smooth'] = spline.derivative()(T)
        
        df.to_csv(output_path, index=False)
        return True, f"Processed with N={N}, added dC1/dT"
    except Exception as e:
        return False, f"Error: {str(e)}"

def process_directory(input_dir, output_dir):
    """Process all CSV files in a directory with automatic N detection."""
    os.makedirs(output_dir, exist_ok=True)
    
    processed_count = 0
    error_count = 0
    skipped_count = 0
    
    print(f"Processing CSV files in: {input_dir}")
    print(f"Saving results to: {output_dir}\n")

    for filename in os.listdir(input_dir):
        if filename.endswith('.csv'):
            input_path = os.path.join(input_dir, filename)
            N = extract_N_from_filename(filename)
            
            if N is None:
                skipped_count += 1
                print(f"[SKIPPED] {filename} - No N value detected in filename")
                continue
                
            output_path = os.path.join(output_dir, filename)
            success, message = process_csv_file(input_path, output_path, N)
            
            if success:
                processed_count += 1
                print(f"[SUCCESS] {filename} - {message}")
            else:
                error_count += 1
                print(f"[ERROR] {filename} - {message}")
    
    print(f"\nProcessing complete!")
    print(f"Successfully processed: {processed_count} files")
    print(f"Skipped (no N value): {skipped_count} files")
    print(f"Errors encountered: {error_count} files")

if __name__ == "__main__":
    INPUT_DIR = "/home/amitkumar/Documents/mpconcmz/build/quadraticpotential"
    OUTPUT_DIR = os.path.join(INPUT_DIR, "processed_csvs")
    
    process_directory(INPUT_DIR, OUTPUT_DIR)