import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import sys
import os

matplotlib.rcParams['font.family'] = 'DejaVu Sans'

if len(sys.argv) != 2:
    print("Usage: python script.py <file_path>")
    sys.exit(1)

file_path = sys.argv[1]
data = pd.read_csv(file_path)

category_counts = data['majority_voting'].value_counts()
category_percentages = category_counts / category_counts.sum() * 100

top_two = category_percentages.sort_values(ascending=False).index.tolist()

top_category = top_two[1] if top_two[0] == "Unassigned" and len(top_two) > 1 else top_two[0]
top_category_clean = str(top_category).replace(" ", "_")

fig_height = 1.5
fig_width = max(6, len(category_percentages) * 2.5)
fig, ax = plt.subplots(figsize=(fig_width, fig_height))

colors = ['#FF9999', '#66B2FF', '#99FF99', '#FFD700', '#FFB6C1', '#8A2BE2']

cumulative = 0
for i, (category, percentage) in enumerate(category_percentages.items()):
    ax.barh([0], [percentage], left=cumulative, color=colors[i % len(colors)],
            label=f"{category}: {percentage:.1f}%")
    cumulative += percentage

ax.axis('off')
ax.legend(loc='upper center',
          bbox_to_anchor=(0.5, 1.2),
          ncol=len(category_percentages),
          frameon=False,
          fontsize=10)
plt.subplots_adjust(left=0.05, right=0.95, top=0.8, bottom=0.2)

file_name_without_extension = os.path.splitext(os.path.basename(file_path))[0]
output_file = f'{file_name_without_extension}_category_percentage_plot_{top_category_clean}.pdf'
plt.savefig(output_file, format='pdf', dpi=300, bbox_inches='tight')
plt.close()

