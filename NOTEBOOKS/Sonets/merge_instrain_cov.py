#!/usr/bin/python
import os
import pandas as pd

# Initialize an empty DataFrame to store merged data
merged_df = pd.DataFrame()

for folder in os.listdir('.'):
    if os.path.isdir(folder) and 'genomeCoverage.csv' in os.listdir(folder):
        # Extract sample name from folder name
        sample_name = folder.split('.')[0]
        
        # Read CSV file (using 1st and 5th columns)
        df = pd.read_csv(
            os.path.join(folder, 'genomeCoverage.csv'),
            usecols=[0, 3],  # Select first and fifth columns
            index_col=0,      # Set first column as index
            header=0         # Use existing header
        )
        
        # Rename fifth column to sample name
        df = df.rename(columns={df.columns[0]: sample_name})
        
        # Merge with main DataFrame
        if merged_df.empty:
            merged_df = df
        else:
            merged_df = pd.merge(merged_df, df, left_index=True, right_index=True, how='outer')

# Save merged DataFrame to CSV
merged_df.to_csv('merged_true_coverage.csv')
