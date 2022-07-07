import pandas as pd
import re
import webbrowser
import time
import random

def clean_csv(file):
    df = pd.read_csv(file)
    df = df.drop(df[df.CDS == 0].index)
    df = df.drop(df[df.tRNA == 0].index)
    df.dropna(subset=['Replicons'], inplace=True)
    df = df.drop_duplicates(subset=['#Organism Name'])  # not sure if I can do it
    df = df.reset_index(drop=True)
    for index, row in df.iterrows():
        try:

            protein_name = re.search(r'chromosome(.*?):(.+?)[\/;]', row['Replicons']).group(2)
            df.loc[index, 'Replicons'] = protein_name

        except AttributeError:
            protein_name = re.search(r'chromosome(.*):(.+)', row['Replicons']).group(2)
            df.loc[index, 'Replicons'] = protein_name

    return df


def create_files_based_gc(pandas_data_frame, num_of_organisms):
    """
    creates 3 files withs NCBI accesion numbers based on gc content
    """
    count_low = 0
    count_medium = 0
    count_high = 0
    low = []
    medium = []
    high = []

    for index, row in pandas_data_frame.iterrows():
        print(row['Replicons'])
        if count_low <= num_of_organisms and row['GC%'] < 40:
            low.append(row['Replicons'])
            count_low += 1

        if count_medium <= num_of_organisms and row['GC%'] in range(40, 60):
            medium.append(row['Replicons'])
            count_medium += 1

        if count_high <= num_of_organisms and row['GC%'] >= 60:
            high.append(row['Replicons'])
            count_high += 1

    for file in ['low', 'medium', 'high']:
        with open(f'to_download_GC/{file}.txt', 'a') as gc:
            if file == 'low':
                sample = random.sample(low, k=100)
            elif file == 'medium':
                sample = random.sample(medium, k=100)
            else:
                sample = random.sample(high, k=100)
            for an in sample:
                gc.write(an + '\n')


def generate_links(file_with_acc_nums):
    with open(file_with_acc_nums) as file:
        for id in file.readlines():
            webbrowser.open(
                f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=<id, {id}, id>6&rettype=fasta_cds_na&retmode=text')
            time.sleep(45)


if __name__ == "__main__":
    files = ["sequences_of_model_organisms/prokaryotes.csv"]
    files_grouped = ['high.txt', 'medium.txt', 'low.txt']

    # create_files_based_gc(clean_csv(files[0]), 1000)

    generate_links(f'to_download_GC/high.txt')
