import sqlite3
import pandas as pd
from parameters import *
import os

def _pre_process_OMA_xref_file(oma_ensembl_file: str, output_processed_file: str,
                                       is_remove_version: bool = True):
    file_tmp = open(output_processed_file, 'a')
    with open(oma_ensembl_file, 'r') as file:
        for line in file:
            if line.startswith("#"):
                if "Format:" in line:
                    new_line = line.split(':')[1].strip().replace("<tab>", '\t')
                else:
                    continue
            elif is_remove_version:
                new_line = line.split('.')[0]
            else:
                new_line = line
            file_tmp.write(new_line.strip() + "\n")
        file.close()
    file_tmp.close()

if __name__ == '__main__':
  sqlite_db = output_file_name + ".db"
  oma_ensembl_processed_file = os.path.join("./", os.path.split(oma_ensembl_file)[1].replace('.txt', '.tsv'))
  oma_entrez_processed_file = os.path.join("./", os.path.split(oma_entrez_file)[1].replace('.txt', '.tsv'))
  _pre_process_OMA_xref_file(oma_ensembl_file, oma_ensembl_processed_file)
  _pre_process_OMA_xref_file(oma_entrez_file, oma_entrez_processed_file, False)
  files_dict = {'oma_entrez': oma_entrez_processed_file,
                'oma_ensembl': oma_ensembl_processed_file,
                'ncbi_ensembl': gene_ensembl_file,
                'bgee_genes': bgee_genes_file}
  if drop_database and os.path.isfile(sqlite_db):
   os.remove(sqlite_db)
  conn = sqlite3.connect(sqlite_db)
  cursor = conn.cursor()
  for table_name, file_path in files_dict.items():
   if file_path.endswith(".csv"):
       data_frame = pd.read_csv(file_path)
   if file_path.endswith(".tsv"):
       data_frame = pd.read_csv(file_path, sep='\t')
   data_frame.to_sql(table_name, conn, index=False, if_exists='replace')
  response_df = pd.read_sql(sql_query, conn)
  output_tsv = open(output_file_name + ".tsv", "w")
  output_tsv.write(response_df.to_csv(sep='\t', index=False))
  output_tsv.close()
  os.remove(oma_ensembl_processed_file)
  os.remove(oma_entrez_processed_file)
  os.remove(sqlite_db)





