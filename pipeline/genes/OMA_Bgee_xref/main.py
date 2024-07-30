import sqlite3
import pandas as pd
from parameters import *
import os
import gzip
import urllib.request


def decompress_gz(file_url: str):
    downloaded_file = os.path.join('./', file_url.rsplit('/', 1)[-1])
    urllib.request.urlretrieve(file_url, downloaded_file)
    if os.path.isfile(downloaded_file):
        with gzip.open(downloaded_file, 'rb') as file:
            file_content = file.readlines()
            os.remove(downloaded_file)
            return file_content
    else:
        raise Exception('The file does not exist: ' + downloaded_file)

def _pre_process_OMA_xref_file(oma_url: str, is_remove_version: bool = True):
    file_output_path = os.path.join('./', oma_url.rsplit('/', 1)[-1].replace('.txt.gz', '.tsv'))
    file_tmp = open(file_output_path, 'a')
    file_content_lines = decompress_gz(oma_url)
    for line in file_content_lines:
        line = line.decode()
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
    file_tmp.close()
    return file_output_path

if __name__ == '__main__':
  sqlite_db = output_file_path + ".db"
  gene_ensembl_file_path = os.path.join('./', 'gene2ensembl.tsv')
  with open(gene_ensembl_file_path, 'wb') as gene_ensembl_file:
      gene_ensembl_file.writelines(decompress_gz(gene_ensembl_url))
  files_dict = {'oma_entrez':  _pre_process_OMA_xref_file(oma_entrez_url, False),
                'oma_ensembl':  _pre_process_OMA_xref_file(oma_ensembl_url),
                'ncbi_ensembl': gene_ensembl_file_path,
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
  if output_file_path.endswith('.tsv'):
      output_tsv = open(output_file_path, "w")
  else:
      output_tsv = open(output_file_path + ".tsv", "w")
  output_tsv.write(response_df.to_csv(sep='\t', index=False))
  output_tsv.close()
  os.remove(files_dict['oma_entrez'])
  os.remove(files_dict['oma_ensembl'])
  os.remove(files_dict['ncbi_ensembl'])
  os.remove(sqlite_db)





