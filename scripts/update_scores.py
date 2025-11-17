import pandas as pd

list_file = 'list'

with open(list_file, 'r') as file:
    new_names = [line.strip() for line in file]

df = pd.read_csv('combined_scores.csv', index_col='name')

unique_new_names = list(set(new_names) - set(df.index))

new_data = pd.DataFrame(index=unique_new_names)
result_df = pd.concat([new_data, df]).fillna('NA')

result_df.to_csv('combined_scores_updated.csv', index_label='name')

print("Updated table saved to combined_scores_updated.csv")
