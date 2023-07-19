## produce test/train in wide format

# ------------- necessary imports
import os
import argparse
import pandas as pd
from sklearn.model_selection import train_test_split


# ------------- code snippet
merged_wide = df.copy()     # creates copy to avoid modifying primary one

for i in range(1, 3):
    train, test = train_test_split (df, test_size=0.2) # split test/train
    train[f'set{i}'] = "train" # assign column names
    test[f'set{i}'] = "test"
        # merges vased on sample/contig columns, validates match 1:1 in rows or drops a column
    merged_wide= merged_wide.merge(pd.concat([train, test]), # concatenates train/test subsets
                                on=["sample", "contig"],
                                validate="1:1")
    
    merged_wide['sample']
    df2 = pd.concat([train, test])
    df2.sort_index(); df2['sample']
