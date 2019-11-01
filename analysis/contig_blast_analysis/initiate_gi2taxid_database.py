import pandas as pd
import s3fs
import multiprocessing as mp
from sqlalchemy import create_engine
engine = create_engine('sqlite:///../../data/gi2taxid.db', echo=False)

gi2taxid_s3_file = "s3://czbiohub-mosquito/blast/gi_to_taxid.txt"

dfs = []
chunksize = 10 ** 6

for chunk in pd.read_csv(gi2taxid_s3_file, header=None, sep=" ", chunksize=chunksize, names=["gi","acc", "taxid"]):
    dfs.append(chunk)
    chunk.to_sql("gi2taxid", con=engine, if_exists="append")
