import pandas as pd
import re
import seaborn as sns
from tabulate import tabulate

fname = "data/result.txt"

sns.set_style('darkgrid')
sns.color_palette('Spectral')

data = pd.read_csv(fname, sep=',', dtype={'rsid':'str', 'chromosome':'object',
                                                           'position':'int', 'genotype':'str'}, comment='#')
df = pd.DataFrame(data)

df['chromosome'] = df['chromosome'].apply(lambda x: re.sub(r'X', r'23', x))
df['chromosome'] = df['chromosome'].apply(lambda x: re.sub(r'MT', r'24', x))
df['chromosome'] = df['chromosome'].apply(lambda x: re.sub(r'Y', r'25', x))
df['chromosome'] = df['chromosome'].apply(lambda x: int(x))

chromosome_dict = {1:'1', 2:'2', 3:'3', 4:'4', 5:'5',
           6:'6', 7:'7', 8:'8', 9:'9', 10:'10',
           11:'11', 12:'12', 13:'13', 14:'14',
           15:'15', 16:'16', 17:'17', 18:'18',
           19:'19', 20:'20', 21:'21', 22:'22',
           23:'X', 24:'MT', 25:'Y'}

df.rename({' rsid': 'rsid'}, axis='columns', inplace=True)

# hiw many SNP's are there per chromosome
rsid_per_chromosome_series = df.groupby('chromosome')['rsid'].count()
rsid_per_chromosome_series.columns = ['chromosome', 'count']
rsid_per_chromosome_series.plot.barh(figsize=(16,9), fontsize=15)
# plt.show()

snp_df = pd.read_csv('data/uniq_snips.csv')

snp_df['genotype'] = snp_df['rsid'].apply(lambda x: re.sub(r'.*([AGCT]);([AGCT])\)', r'\1\2', x))
new_cols = ['rsid', 'magnitude', 'repute', 'summary', 'genotype']
snp_df.columns = new_cols

snp_df['rsid'] = snp_df['rsid'].map(lambda x : x.lower())
snp_df['rsid'] = snp_df['rsid'].map(lambda x : re.sub(r'([a-z]{1,}[\d]+)\([agct];[agct]\)', r'\1', x))

null_repute = snp_df[snp_df['repute'].isnull()]
null_summaries = snp_df[snp_df['summary'].isnull()]
null_repute_and_summaries = pd.concat([null_repute,null_summaries]).drop_duplicates().reset_index(drop=True)

snp_df['repute'].fillna(value='Neutral', inplace=True)
snp_df['summary'].fillna(value='None', inplace=True)

new_df = snp_df.merge(df, how='inner', on=['rsid', 'genotype'], suffixes=('_SNPedia', '_myDNA'))

good_genes = new_df[new_df.repute == 'Good']
bad_genes = new_df[new_df.repute == 'Bad']
interesting_genes = new_df[new_df.magnitude > 4] # 4 is the threshold for "worth your time" given by SNPedia

base_url = 'https://www.snpedia.com/index.php/'

gene_urls = [base_url + rsid for rsid in bad_genes['rsid']]
for url in gene_urls:
    print(url, '\n')

print(tabulate(bad_genes))
