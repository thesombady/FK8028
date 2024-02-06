import pandas as pd


"""
for i in range(0, 10001):

    df = pd.read_csv('iteration/lattice{}.csv'.format(i))
"""
    
df = pd.read_csv('iteration/lattice0.csv')
for i in range(1, 10001):
    df = df.add(pd.read_csv('iteration/lattice{}.csv'.format(i)), fill_value=0)
print(df.loc[2])
