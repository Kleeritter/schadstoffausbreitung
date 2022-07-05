import pandas as pd

df = pd.DataFrame.from_dict({
    'Name': ['Ray', 'Jane', 'Kate', 'Nik', 'Autumn', 'Kasi', 'Mandeep', 'Evan', 'Kyra', 'Jim'],
    'Age': [12, 7, 33, 34, 45, 65, 77, 11, 32, 55]
})

df['Age Groups'] = pd.qcut(df['Age'], 4)
print(df.head())