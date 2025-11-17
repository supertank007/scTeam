import os
import pandas as pd

input_dir = 'cellR'


data = {}

for filename in os.listdir(input_dir):
    filepath = os.path.join(input_dir, filename)
    if os.path.isfile(filepath):
        with open(filepath, 'r') as file:
            lines = file.readlines()
            for line in lines:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    name = parts[0]
                    score = parts[1]
                    if filename not in data:
                        data[filename] = {}
                    data[filename][name] = score

result_df = pd.DataFrame(data).T.fillna('NA')

result_df.to_csv('combined_scores.csv', index_label='name')
print("Combined table saved to combined_scores.csv")
