from human_model import process_dna, model
from argparse import ArgumentParser
import os

argparse = ArgumentParser()
argparse.add_argument(
    "-i",
    "--infile",
    help="Path to input file with Homo sapiens sequences",
    required=True,
)
args = argparse.parse_args()

inf = args.infile


def predict_genes(infile, model=model):
    X, proteins = process_dna(infile)
    headers = list(X.keys())
    predictions = []
    for x in list(X.values()):
        p = model.predict(x)
        predictions.append(p)
        dirn = os.path.dirname(infile)
        bsnm = os.path.splitext(os.path.basename(infile))[0]
        new_file = os.path.join(dirn, bsnm + ".genecalled.csv")
        file_obj = open(new_file, "w")
        file_obj.write(f"ORF_ID,PROTEIN_SEQ,PREDICTED_BIOTYPE\n")
    for i in range(len(predictions)):
        file_obj.write(f"{headers[i]},{proteins[headers[i]]},{predictions[i][0]}\n")


if __name__ == "__main__":
    predict_genes(inf)
