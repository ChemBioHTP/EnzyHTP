"""This module contains helper functions related to Pandas library (https://pandas.pydata.org/)
Pandas is widely used in EnzyHTP for PDB parsing etc.

Author: QZ Shao, <shaoqz@icloud.com>
Date: 2022-10-21
"""
from typing import List
import pandas as pd

def split_df_base_on_column_value(df: pd.DataFrame,
                                  column_name: str,
                                  split_values: list,
                                  copy: bool = False) -> List[pd.DataFrame]:
    """
    split a dataframe base on the value of a column
    ** the line in the split values will not be included **
    Arg:
        df: the target dataframe
        column_name: the reference column"s name
        split_value: the value mark for spliting
    """
    # empty list
    if not split_values:
        if copy:
            return [df.copy()]
        return [df]

    split_values = sorted(split_values)
    frist = 1
    result_dfs = []
    for this_value in split_values:
        if frist:
            frist_df = df[df[column_name] < this_value]
            result_dfs.append(frist_df)
            frist = 0
            last_value = this_value
            continue
        result_df = df[(df[column_name] < this_value) & (df[column_name] > last_value)]
        result_dfs.append(result_df)
        last_value = this_value
    # deal with the last portion
    result_dfs.append(df[df[column_name] > split_values[-1]])

    if copy:
        for i, df_i in enumerate(result_dfs):
            result_dfs[i] = df_i.copy()

    return result_dfs


def batch_edit_df_loc_value(df: pd.DataFrame, loc_value_list: List[tuple], column: str):
    """
    batch edit "column" of "df" with the "loc_value_list"
    """
    for loc, value in loc_value_list:
        df.loc[loc, column] = value
