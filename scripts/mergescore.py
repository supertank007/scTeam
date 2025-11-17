import pandas as pd

df1 = pd.read_csv('file1.csv', index_col=0)
df2 = pd.read_csv('file2.csv', index_col=0)

# 打印两个 DataFrame 的索引（行名）
print("df1 index:", df1.index)
print("df2 index:", df2.index)

# 按行名（索引）合并两个表格
merged_df = df1.join(df2, how='outer')  # 可以根据需要使用 'inner', 'left', 'right' 或 'outer'

# 打印合并后的表格
print(merged_df)

# 如果你希望将合并后的表格保存为新的 CSV 文件
merged_df.to_csv('merged_file.csv')

import pandas as pd

# 读取 CSV 文件
df = pd.read_csv('your_file.csv')

subset = df.iloc[:, 9:]

# 删除这些列全是空值的行
df = df[~subset.isnull().all(axis=1)]

# 打印结果
print(df)

# 如果你希望将处理后的表格保存为新的 CSV 文件
df.to_csv('cleaned_file.csv', index=False)


import pandas as pd

# 读取处理后的 CSV 文件
df = pd.read_csv('cleaned_file.csv')

# 获取第7列（分组信息）和第10列开始的所有列（数据）
group_info = df.iloc[:, 6]
data = df.iloc[:, 9:]

# 合并分组信息和数据
df = pd.concat([group_info, data], axis=1)

# 分组计算每列的求和
grouped_sum = df.groupby(df.columns[0]).sum()

# 计算所有分组的列总和
total_sum = grouped_sum.sum(axis=0)

# 计算每个分组的比例
proportions = grouped_sum.div(total_sum, axis=1)

# 打印结果
print(proportions)

# 保存计算后的比例
proportions.to_csv('group_proportions.csv')


